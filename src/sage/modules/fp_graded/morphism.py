# sage.doctest: needs sage.combinat
r"""
Homomorphisms of finitely presented graded modules

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

# *****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu>
#             and          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from inspect import isfunction

from sage.categories.homset import End, Hom
from sage.categories.morphism import Morphism
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.infinity import infinity


def _create_relations_matrix(module, relations, source_degs, target_degs):
    r"""
    The action by the given relations can be written as multiplication by
    the matrix `R = (r_{ij})_{i,j}` where each entry is an algebra element and
    each row in the matrix contains the coefficients of a single relation.

    For a given source degree, `n`, the multiplication by `r_{ij}` restricts to
    a linear transformation `M_n\to M_{n + \deg(r_{ij})}`.  This function returns
    the matrix of linear transformations gotten by restricting `R` to the given
    source degrees.

    INPUT:

    - ``module`` -- the module where the relations acts
    - ``relations`` -- list of lists of algebra coefficients defining the
      matrix `R`
    - ``source_degs`` -- list of integer degrees; its length should be
      equal to the number of columns of `R`
    - ``target_degs`` -- list of integer degrees; its length should be
      equal to the number of rows of `R`

    Furthermore must the degrees given by the input satisfy the following::

        source_degs[j] + deg(r[i,j]) = target_degs[i]

    for all `i, j`.

    OUTPUT:

    - ``block_matrix`` -- list of lists representing a matrix of linear
      transformations `(T_{ij})`.  Each transformtion `T_{ij}` is the linear map
      representing multiplication by the coefficient `r_{ij}` restricted to
      the module elements of degree ``source_degs[j]``.
    - ``R`` -- a matrix representing ``block_matrix`` as a single linear
      transformation

    TESTS::

        sage: from sage.modules.fp_graded.module import FPModule
        sage: from sage.modules.fp_graded.morphism import _create_relations_matrix
        sage: A = SteenrodAlgebra(p=2)
        sage: M = FPModule(A, [0], [[0]])
        sage: blocks, R = _create_relations_matrix(M, [[Sq(2)]], [4], [6])

        sage: blocks
          [[Vector space morphism represented by the matrix:
            [0 1 0]
            [0 1 1]
            Domain: Vector space quotient V/W of dimension 2 over Finite Field of size 2 where
            V: Vector space of dimension 2 over Finite Field of size 2
            W: Vector space of degree 2 and dimension 0 over Finite Field of size 2
            Basis matrix:
            []
            Codomain: Vector space quotient V/W of dimension 3 over Finite Field of size 2 where
            V: Vector space of dimension 3 over Finite Field of size 2
            W: Vector space of degree 3 and dimension 0 over Finite Field of size 2
            Basis matrix:
            []]]

        sage: R
          [0 0]
          [1 1]
          [0 1]
    """
    from sage.matrix.constructor import matrix

    if not relations:
        raise ValueError('no relations given, can not build matrix')

    # Create the block matrix of linear transformations.
    block_matrix = []
    for i, r_i in enumerate(relations):
        row = []
        target_space = module.vector_presentation(target_degs[i])

        for j, r_ij in enumerate(r_i):

            values = []
            for b in module.basis_elements(source_degs[j]):
                w = r_ij * b
                values.append(
                    target_space.zero() if w.is_zero() else w.vector_presentation())

            row.append(
                Hom(module.vector_presentation(source_degs[j]), target_space)(values))

        block_matrix.append(row)

    # Deal with the case of zero dimensional matrices first.
    total_source_dim = 0
    for el in block_matrix[0]:
        total_source_dim += el.domain().dimension()
    total_target_dim = 0
    for row in block_matrix:
        total_target_dim += row[0].codomain().dimension()
    if total_source_dim == 0:
        return block_matrix, matrix(total_target_dim, 0)
    elif total_target_dim == 0:
        return block_matrix, matrix(0, total_source_dim)

    # Create a matrix from the matrix of linear transformations.
    entries = []
    for j, block in enumerate(block_matrix[0]):
        for n in range(block.domain().dimension()):
            column = []
            for i in range(len(block_matrix)):
                lin_trans = block_matrix[i][j]
                column += lin_trans(lin_trans.domain().basis()[n])
            entries.append(column)

    return block_matrix, matrix(module.base_ring().base_ring(), entries).transpose()


class FPModuleMorphism(Morphism):
    r"""
    Create a homomorphism between finitely presented graded modules.

    INPUT:

    - ``parent`` -- a homspace of finitely presented graded modules
    - ``values`` -- list of elements in the codomain; each element
      corresponds to a module generator in the domain
    - ``check`` -- boolean (default: ``True``); if ``True``, check
      that the morphism is well-defined

    TESTS::

        sage: from sage.modules.fp_graded.module import FPModule

    Trying to map the generators of a non-free module into a
    free module::

        sage: A = SteenrodAlgebra(2)
        sage: F = FPModule(A, [2,3])
        sage: Q = FPModule(A, [2,3], relations=[[Sq(6), Sq(5)]])
        sage: v1 = Q((Sq(1), 0))
        sage: v2 = Q((0, 1))
        sage: m = Hom(F, Q)( (v1, v2) )
        Traceback (most recent call last):
        ...
        ValueError: ill-defined homomorphism: degrees do not match

    Trying to map the generators of a non-free module into a
    free module::

        sage: w = Hom(Q, F)( (F((1, 0)), F((0, 1))) )
        Traceback (most recent call last):
        ...
        ValueError: relation Sq(6)*g[2] + Sq(5)*g[3] is not sent to zero
    """

    def __init__(self, parent, values, check=True):
        r"""
        Create a homomorphism between finitely presented graded modules.

        OUTPUT:

        A module homomorphism defined by sending the generator
        with index `i` to the corresponding element in ``values``.

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: M = FPModule(A2, [0], relations=[[Sq(1)]])
            sage: N = FPModule(A2, [0], relations=[[Sq(4)],[Sq(1)]])
            sage: f = Hom(M,N)([A2.Sq(3)*N.generator(0)])
            sage: TestSuite(f).run()
        """
        from .homspace import FPModuleHomspace

        if not isinstance(parent, FPModuleHomspace):
            raise TypeError('parent (=%s) must be a fp module hom space' % parent)

        # Get the values.
        C = parent.codomain()
        D = parent.domain()
        if isfunction(values):
            values = [C(values(g)) for g in D.generators()]
        elif not values:
            values = len(D.generator_degrees()) * [C.zero()]
        else:
            values = [C(a) for a in values]

        # Check the homomorphism is well defined.
        if len(D.generator_degrees()) != len(values):
            raise ValueError('the number of values must equal the number of '
                             'generators in the domain; invalid argument: %s' % values)

        self._values = tuple(values)

        # Call the base class constructor.
        Morphism.__init__(self, parent)

        # Check that the homomorphism is well defined.
        if check:
            self._free_morphism
            for relation in parent.domain().relations():
                # The relation is an element in the free part of the domain.
                img = self._free_morphism(relation)
                if parent.codomain()(img):
                    raise ValueError('relation %s is not sent to zero' % relation)

    @lazy_attribute
    def _free_morphism(self):
        """
        Return the lifted morphism between free modules.

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: M = FPModule(A2, [0], relations=[[Sq(1)]])
            sage: N = FPModule(A2, [0], relations=[[Sq(4)],[Sq(1)]])
            sage: f = Hom(M,N)([A2.Sq(3)*N.generator(0)])
            sage: f._free_morphism
            Module endomorphism of Free graded left module on 1 generator
             over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[0] |--> Sq(3)*g[0]
        """
        P = self.parent()
        from sage.modules.fp_graded.free_module import FreeGradedModule
        if isinstance(P.codomain(), FreeGradedModule):
            Homspace = Hom(P.domain()._j.codomain(), P.codomain())
            return Homspace(self._values)
        Homspace = Hom(P.domain()._j.codomain(), P.codomain()._j.codomain())
        return Homspace([v.lift_to_free() for v in self._values])

    def change_ring(self, algebra):
        r"""
        Change the base ring of ``self``.

        INPUT:

        - ``algebra`` -- a graded algebra

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: M = FPModule(A2, [0], relations=[[Sq(1)]])
            sage: N = FPModule(A2, [0], relations=[[Sq(4)],[Sq(1)]])

            sage: f = Hom(M,N)([A2.Sq(3)*N.generator(0)]); f
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              To:   Finitely presented left module on 1 generator and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[0] |--> Sq(3)*g[0]

            sage: f.base_ring()
            sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]

            sage: g = f.change_ring(A3)
            sage: g.base_ring()
            sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [4, 3, 2, 1]
        """
        new_codomain = self.codomain().change_ring(algebra)
        # We have to change the ring for the values, too:
        new_values = [new_codomain([algebra(a)
                                    for a in v.dense_coefficient_list()])
                      for v in self._values]
        return Hom(self.domain().change_ring(algebra), new_codomain)(new_values)

    def degree(self):
        r"""
        The degree of ``self``.

        OUTPUT: the integer degree of ``self``

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])
            sage: homspace = Hom(M, N)

            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f.degree()
            7

        The trivial homomorphism has no degree::

            sage: homspace.zero().degree()
            Traceback (most recent call last):
            ...
            ValueError: the zero morphism does not have a well-defined degree

        TESTS::

            sage: M = FPModule(SteenrodAlgebra(p=2), [7])
            sage: N = FPModule(SteenrodAlgebra(p=2), [0], relations=[[Sq(1)]])
            sage: f = Hom(M,N)([Sq(1)*N.generator(0)])
            sage: f == Hom(M,N).zero()
            True
            sage: f.degree()
            Traceback (most recent call last):
            ...
            ValueError: the zero morphism does not have a well-defined degree
        """
        if self.is_zero():
            # The zero morphism has no degree.
            raise ValueError("the zero morphism does not have a well-defined degree")
        return self._free_morphism.degree()

    def values(self):
        r"""
        The values under ``self`` of the module generators of the domain module.

        OUTPUT: a sequence of module elements of the codomain

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])
            sage: homspace = Hom(M, N)

            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)

            sage: f.values()
            (Sq(5)*g[2], Sq(3,1)*g[2])

            sage: homspace.zero().values()
            (0, 0)

            sage: homspace = Hom(A.free_graded_module((0,1)), A.free_graded_module((2,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f.values()
            (Sq(5)*g[2], Sq(3,1)*g[2])
            sage: homspace.zero().values()
            (0, 0)
        """
        return self._values

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])
            sage: homspace = Hom(M, N)
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f._richcmp_(f, op=2)
            True
            sage: f._richcmp_(f, op=3)
            False

            sage: homspace = Hom(A.free_graded_module((0,1)), A.free_graded_module((2,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f._richcmp_(f, op=2)
            True
            sage: f._richcmp_(f, op=3)
            False
        """
        try:
            same = (self - other).is_zero()
        except ValueError:
            return False

        # Equality
        if op == 2:
            return same

        # Non-equality
        if op == 3:
            return not same

        return False

    def __add__(self, g):
        r"""
        The pointwise sum of ``self`` and ``g``.

        Pointwise addition of two homomorphisms `f` and `g` with the same domain
        and codomain is given by the formula `(f+g)(x) = f(x) + g(x)` for
        every `x` in the domain of `f`.

        INPUT:

        - ``g`` -- a homomorphism with the same domain and codomain as ``self``

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])
            sage: v = N.generator(0)
            sage: homspace = Hom(M, N)
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = homspace(values)
            sage: ff = f + f
            sage: ff.is_zero()
            True
            sage: ff + f == f
            True

            sage: N = A.free_graded_module((2,))
            sage: v = N.generator(0)
            sage: homspace = Hom(A.free_graded_module((0,1)), N)
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = homspace(values)
            sage: ff = f + f
            sage: ff.is_zero()
            True
            sage: ff + f == f
            True
            sage: ff = f + f
            sage: ff.is_zero()
            True
        """
        if self.domain() != g.domain():
            raise ValueError('morphisms do not have the same domain')
        elif self.codomain() != g.codomain():
            raise ValueError('morphisms do not have the same codomain')
        if self.is_zero():
            return g
        if g.is_zero():
            return self
        if self.degree() and g.degree() and self.degree() != g.degree():
            raise ValueError('morphisms do not have the same degree')

        v = [self(x) + g(x) for x in self.domain().generators()]

        return self.parent()(v)

    def __neg__(self):
        r"""
        The additive inverse of ``self`` with respect to the group
        structure given by pointwise sum.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])
            sage: v = N.generator(0)
            sage: homspace = Hom(M, N)
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = homspace(values)
            sage: f_neg = -f; f_neg
            Module morphism:
              From: Finitely presented left module on 2 generators and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(5)*g[2]
                    g[1] |--> Sq(3,1)*g[2]
            sage: (f + f_neg).is_zero()
            True

            sage: N = A.free_graded_module((2,))
            sage: v = N.generator(0)
            sage: homspace = Hom(A.free_graded_module((0,1)), N)
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = homspace(values)
            sage: f_neg = -f; f_neg
            Module morphism:
              From: Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(5)*g[2]
                    g[1] |--> Sq(3,1)*g[2]
            sage: (f + f_neg).is_zero()
            True
        """
        return self.parent()([-x for x in self._values])

    def __sub__(self, g):
        r"""
        The difference between ``self`` and ``g`` with
        respect to the group structure given by pointwise sum.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0])
            sage: N = FPModule(A, [0], [[Sq(4)]])
            sage: f = Hom(M, N)( [Sq(3)*N.generator(0)] )
            sage: g = Hom(M, N)( [Sq(0,1)*N.generator(0)] )
            sage: f - g
            Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> (Sq(0,1)+Sq(3))*g[0]

            sage: f = Hom(M, N)( [Sq(4)*N.generator(0)] ) # the zero map
            sage: g = Hom(M, N)( [Sq(1,1)*N.generator(0)] )
            sage: f - g
            Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(1,1)*g[0]

            sage: N = A.free_graded_module((2,))
            sage: v = N.generator(0)
            sage: homspace = Hom(A.free_graded_module((0,1)), N)
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = homspace(values)
            sage: values2 = [Sq(5)*v, Sq(3,1)*v]
            sage: g = homspace(values2)
            sage: f - g == 0
            True
        """
        return self.__add__(g.__neg__())

    def __mul__(self, g):
        r"""
        The composition of ``g`` followed by ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0], [[Sq(1,2)]])
            sage: N = FPModule(A, [0], [[Sq(2,2)]])
            sage: f = Hom(M, N)( [Sq(2)*N.generator(0)] )
            sage: g = Hom(N, M)( [Sq(2,2)*M.generator(0)] )
            sage: fg = f * g; fg
            Module endomorphism of Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> (Sq(0,1,1)+Sq(1,3)+Sq(3,0,1))*g[0]
            sage: fg.is_endomorphism()
            True

            sage: M = A.free_graded_module((0,1))
            sage: N = A.free_graded_module((2,))
            sage: v = N.generator(0)
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = Hom(M, N)(values)
            sage: values2 = [Sq(2)*M.generator(0)]
            sage: g = Hom(N, M)(values2)
            sage: fg = f * g; fg
            Module endomorphism of Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[2] |--> (Sq(4,1)+Sq(7))*g[2]
            sage: fg.is_endomorphism()
            True

        TESTS::

            sage: from sage.modules.fp_graded.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, (0,1))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f * f
            Traceback (most recent call last):
            ...
            ValueError: morphisms not composable

            sage: M = A.free_graded_module((0,1))
            sage: N = A.free_graded_module((2,))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f * f
            Traceback (most recent call last):
            ...
            ValueError: morphisms not composable
        """
        if self.parent().domain() != g.parent().codomain():
            raise ValueError('morphisms not composable')
        homset = Hom(g.parent().domain(), self.parent().codomain())
        return homset([self(g(x)) for x in g.domain().generators()])

    @cached_method
    def is_zero(self):
        r"""
        Decide if ``self`` is the zero homomorphism.

        OUTPUT: the boolean value ``True`` if ``self`` is trivial and ``False`` otherwise

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]

            sage: f = Hom(M, N)(values)
            sage: f.is_zero()
            False

            sage: (f-f).is_zero()
            True

            sage: M = A.free_graded_module((0,1))
            sage: N = A.free_graded_module((2,))
            sage: v = N.generator(0)
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = Hom(M, N)(values)
            sage: f.is_zero()
            False
            sage: (f-f).is_zero()
            True
        """
        return all(x.is_zero() for x in self._values)

    __bool__ = is_zero

    @cached_method
    def is_identity(self):
        r"""
        Decide if ``self`` is the identity endomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]

            sage: f = Hom(M, N)(values)
            sage: f.is_identity()
            False

            sage: one = Hom(M, M)(M.generators()); one
            Module endomorphism of Finitely presented left module on 2 generators and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> g[0]
                    g[1] |--> g[1]

            sage: one.is_identity()
            True

            sage: M = A.free_graded_module((0,1))
            sage: N = A.free_graded_module((2,))
            sage: v = N.generator(0)
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = Hom(M, N)(values)
            sage: f.is_identity()
            False
            sage: one = Hom(M, M)(M.generators()); one
            Module endomorphism of Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> g[0]
                    g[1] |--> g[1]
            sage: one.is_identity()
            True
        """
        return (self.parent().is_endomorphism_set()
                and self.parent().identity() == self)

    def __call__(self, x):
        r"""
        Evaluate the homomorphism at the given domain element ``x``.

        INPUT:

        - ``x`` -- an element of the domain of the homomorphism

        OUTPUT:

        The module element of the codomain which is the value of ``x``
        under this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])

            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)

            sage: f.__call__(M.generator(0))
            Sq(5)*g[2]

            sage: f.__call__(M.generator(1))
            Sq(3,1)*g[2]
        """
        if x.parent() != self.domain():
            raise ValueError('cannot evaluate morphism on element not in domain')

        return self.codomain()(self._free_morphism(x.lift_to_free()))

    def _repr_type(self):
        """
        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: M = FPModule(A2, [0], relations=[[Sq(1)]])
            sage: N = FPModule(A2, [0], relations=[[Sq(4)],[Sq(1)]])
            sage: f = Hom(M,N)([A2.Sq(3)*N.generator(0)])
            sage: f._repr_type()
            'Module'
        """
        return "Module"

    def _repr_defn(self):
        """
        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: M = FPModule(A2, [0], relations=[[Sq(1)]])
            sage: N = FPModule(A2, [0], relations=[[Sq(4)],[Sq(1)]])
            sage: f = Hom(M,N)([A2.Sq(3)*N.generator(0)])
            sage: f._repr_defn()
            'g[0] |--> Sq(3)*g[0]'

            sage: A = SteenrodAlgebra(2)
            sage: F1 = A.free_graded_module((4,5), names='b')
            sage: F2 = A.free_graded_module((3,4), names='c')
            sage: H = Hom(F1, F2)
            sage: f = H((F2((Sq(4), 0)), F2((0, Sq(4)))))
            sage: print(f._repr_defn())
            b[4] |--> Sq(4)*c[3]
            b[5] |--> Sq(4)*c[4]
        """
        s = '\n'.join(['%s |--> %s' % (x, y) for (x, y) in
                       zip(self.domain().generators(), self._values)])
        return s

    @cached_method
    def vector_presentation(self, n):
        r"""
        Return the restriction of ``self`` to the domain module elements
        of degree ``n``.

        The restriction of a nonzero module homomorphism to the free module
        of module elements of degree `n` is a linear function into the free
        module of elements of degree `n+d` belonging to the codomain.
        Here `d` is the degree of this homomorphism.

        When this homomorphism is zero, it has no well defined degree so the
        function cannot be presented since we do not know the degree of its
        codomain.  In this case, an error is raised.

        INPUT:

        - ``n`` -- integer degree

        OUTPUT:

        A linear function of finite dimensional free modules over the
        ground ring of the algebra for this module.  The domain is isomorphic
        to the free module of domain elements of degree ``n`` of ``self``
        via the choice of basis given by
        :meth:`~sage.modules.fp_graded.free_module.FreeGradedModule.basis_elements`.
        If the morphism is zero, an error is raised.

        .. SEEALSO::

            * :meth:`sage.modules.fp_graded.module.FPModule.vector_presentation`
            * :meth:`sage.modules.fp_graded.module.FPModule.basis_elements`
            * :meth:`sage.modules.fp_graded.free_module.FreeGradedModule.vector_presentation`
            * :meth:`sage.modules.fp_graded.free_module.FreeGradedModule.basis_elements`

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = FPModule(A, [2], [[Sq(4)]])
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.vector_presentation(0)
            Vector space morphism represented by the matrix:
            [0]
            Domain: Vector space quotient V/W of dimension 1 over Finite Field of size 2 where
            V: Vector space of dimension 1 over Finite Field of size 2
            W: Vector space of degree 1 and dimension 0 over Finite Field of size 2
            Basis matrix:
            []
            Codomain: Vector space quotient V/W of dimension 1 over Finite Field of size 2 where
            V: Vector space of dimension 2 over Finite Field of size 2
            W: Vector space of degree 2 and dimension 1 over Finite Field of size 2
            Basis matrix:
            [0 1]
            sage: f.vector_presentation(1)
            Vector space morphism represented by the matrix:
            [0 0]
            [0 1]
            Domain: Vector space quotient V/W of dimension 2 over Finite Field of size 2 where
            V: Vector space of dimension 2 over Finite Field of size 2
            W: Vector space of degree 2 and dimension 0 over Finite Field of size 2
            Basis matrix:
            []
            Codomain: Vector space quotient V/W of dimension 2 over Finite Field of size 2 where
            V: Vector space of dimension 3 over Finite Field of size 2
            W: Vector space of degree 3 and dimension 1 over Finite Field of size 2
            Basis matrix:
            [0 1 1]
            sage: f.vector_presentation(2)
            Vector space morphism represented by the matrix:
            [0 0]
            Domain: Vector space quotient V/W of dimension 1 over Finite Field of size 2 where
            V: Vector space of dimension 2 over Finite Field of size 2
            W: Vector space of degree 2 and dimension 1 over Finite Field of size 2
            Basis matrix:
            [1 1]
            Codomain: Vector space quotient V/W of dimension 2 over Finite Field of size 2 where
            V: Vector space of dimension 4 over Finite Field of size 2
            W: Vector space of degree 4 and dimension 2 over Finite Field of size 2
            Basis matrix:
            [0 0 1 0]
            [0 0 0 1]

            sage: M = A.free_graded_module((0,1))
            sage: N = A.free_graded_module((2,))
            sage: v = N.generator(0)
            sage: values = [Sq(5)*v, Sq(3,1)*v]
            sage: f = Hom(M, N)(values)
            sage: f.vector_presentation(0)
            Vector space morphism represented by the matrix:
            [0 1]
            Domain: Vector space of dimension 1 over Finite Field of size 2
            Codomain: Vector space of dimension 2 over Finite Field of size 2
            sage: f.vector_presentation(1)
            Vector space morphism represented by the matrix:
            [0 0 0]
            [0 1 0]
            Domain: Vector space of dimension 2 over Finite Field of size 2
            Codomain: Vector space of dimension 3 over Finite Field of size 2
            sage: f.vector_presentation(2)
            Vector space morphism represented by the matrix:
            [0 0 1 1]
            [0 0 0 0]
            Domain: Vector space of dimension 2 over Finite Field of size 2
            Codomain: Vector space of dimension 4 over Finite Field of size 2

        TESTS::

            sage: F = FPModule(A, [0])
            sage: Q = FPModule(A, [0], [[Sq(2)]])
            sage: z = Hom(F, Q)([Sq(2)*Q.generator(0)])
            sage: z.is_zero()
            True
            sage: z.vector_presentation(0)
            Traceback (most recent call last):
            ...
            ValueError: the zero map has no vector presentation
        """
        # The trivial map has no degree, so we can not create the codomain
        # of the linear transformation.
        if self.is_zero():
            raise ValueError("the zero map has no vector presentation")

        D_n = self.domain().vector_presentation(n)
        C_n = self.codomain().vector_presentation(self.degree() + n)

        values = [self(e) for e in self.domain().basis_elements(n)]

        return Hom(D_n, C_n)([C_n.zero() if e.is_zero() else e.vector_presentation()
                              for e in values])

    def solve(self, x):
        r"""
        Return an element in the inverse image of ``x``.

        INPUT:

        - ``x`` -- an element of the codomain of this morphism

        OUTPUT:

        An element of the domain which maps to ``x`` under this morphism
        or ``None`` if ``x`` was not in the image of this morphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0], [[Sq(3)]])
            sage: N = FPModule(A, [0], [[Sq(2,2)]])
            sage: f = Hom(M, N)( [Sq(2)*N.generator(0)] )
            sage: y = Sq(1,1)*N.generator(0); y
            Sq(1,1)*g[0]
            sage: x = f.solve(y); x
            Sq(2)*g[0]
            sage: y == f(x)
            True

        Trying to lift an element which is not in the image results
        in a ``None`` value::

            sage: z = f.solve(Sq(1)*N.generator(0))
            sage: z is None
            True

        TESTS::

            sage: f.solve(Sq(2,2)*M.generator(0))
            Traceback (most recent call last):
            ...
            ValueError: the given element is not in the codomain of this homomorphism
        """
        if x.parent() != self.codomain():
            raise ValueError('the given element is not in the codomain of this homomorphism')

        # The zero element lifts over all morphisms.
        if x.is_zero():
            return self.domain().zero()

        # Handle the trivial homomorphism since it does not have a well defined
        # degree.
        if self.is_zero():
            return None

        # Handle the case where both morphism and element are non-trivial.
        n = x.degree() - self.degree()
        f_n = self.vector_presentation(n)

        v = x.vector_presentation()

        # Return None if ``x`` cannot be lifted.
        if v not in f_n.image():
            return None

        u = f_n.matrix().solve_left(v)
        return self.domain().element_from_coordinates(u, n)

    def lift(self, f, verbose=False):
        r"""
        Return a lift of this homomorphism over the given homomorphism ``f``.

        INPUT:

        - ``f`` -- a homomorphism with codomain equal to the codomain of ``self``
        - ``verbose`` -- boolean (default: ``False``); enable progress messages

        OUTPUT:

        A homomorphism `g` with the property that ``self`` equals `f \circ g`.
        If no lift exist, ``None`` is returned.

        ALGORITHM:

        Let `s` be this homomorphism with `L` the domain of `s`.
        Choose `x_1, \ldots, x_N` such that `f(x_i) = s(g_i)`,
        where the `g_i`'s are the module generators of `L`.

        The linear function sending `g_i` to `x_i` for every `i` is well
        defined if and only if the vector `x = (x_1, \ldots, x_N)` lies
        in the nullspace of the coefficient matrix `R = (r_{ij})` given
        by the relations of `L`.

        Let `k \in \ker(f)` solve the matrix equation:

        .. MATH::

            R \cdot k = R \cdot x.

        Define a module homomorphism by sending the generators of `L` to
        `x_1 - k_1, \ldots, x_N - k_N`.  This is well defined, and is also a
        lift of this homomorphism over `f`.

        Note that it does not matter how we choose the initial elements `x_i`:
        If `x'` is another choice then `x' - x\in \ker(f)` and
        `R\cdot k = R\cdot x` if and only if `R\cdot (k + x' - x) = R\cdot x'`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)

        Lifting a map from a free module is always possible::

            sage: M = FPModule(A, [0], [[Sq(3)]])
            sage: N = FPModule(A, [0], [[Sq(2,2)]])
            sage: F = FPModule(A, [0])
            sage: f = Hom(M,N)([Sq(2)*N.generator(0)])
            sage: k = Hom(F,N)([Sq(1)*Sq(2)*N.generator(0)])
            sage: f_ = k.lift(f)
            sage: f*f_ == k
            True
            sage: f_
            Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(1)*g[0]

        A split projection::

            sage: A_plus_HZ = FPModule(A, [0,0], [[0, Sq(1)]])
            sage: HZ = FPModule(A, [0], [[Sq(1)]])
            sage: q = Hom(A_plus_HZ, HZ)([HZ([1]), HZ([1])])
            sage: # We can construct a splitting of `q` manually:
            sage: split = Hom(HZ,A_plus_HZ)([A_plus_HZ.generator(1)])
            sage: (q*split).is_identity()
            True

        Thus, lifting the identity homomorphism over `q` should be possible::

            sage: one = Hom(HZ,HZ).identity()
            sage: j = one.lift(q); j
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> g[0, 1]
            sage: q*j
            Module endomorphism of Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> g[0]

        Lifting over the inclusion of the image sub module::

            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0], relations=[[Sq(0,1)]])
            sage: f = Hom(M,M)([Sq(2)*M.generator(0)])
            sage: im = f.image(top_dim=10)
            sage: f.lift(im)
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> g[2]

        When a lift cannot be found, the ``None`` value is returned.  By
        setting the verbose argument to ``True``, an explanation of why
        the lifting failed will be displayed::

            sage: F2 = FPModule(A, [0,0])
            sage: non_surjection = Hom(F2, F2)([F2([1, 0]), F2([0, 0])])
            sage: lift = Hom(F2, F2).identity().lift(non_surjection, verbose=True)
            The generators of the domain of this homomorphism do not map into
             the image of the homomorphism we are lifting over.
            sage: lift is None
            True

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: # The trivial map often involved in corner cases..
            sage: trivial_map = Hom(FPModule(A, [0]), FPModule(A, [])).zero()
            sage: trivial_map.lift(trivial_map) == 0
            True

            sage: F = FPModule(A, [0])
            sage: HZ = FPModule(A, [0], relations=[[Sq(1)]])
            sage: f = Hom(F,HZ)(HZ.generators())
            sage: split = Hom(HZ, HZ).identity().lift(f, verbose=True)
            The homomorphism cannot be lifted in any way such that the relations
             of the domain are respected: matrix equation has no solutions
            sage: split is None
            True

            sage: Hom(F, F).identity().lift(f, verbose=true)
            Traceback (most recent call last):
            ...
            ValueError: the codomains of this homomorphism and the homomorphism
             we are lifting over are different

            sage: f.lift(Hom(HZ, HZ).zero(), verbose=True)
            This homomorphism cannot lift over a trivial homomorphism since
             it is non-trivial.

            sage: Ap = SteenrodAlgebra(p=2, profile=(2,2,2,1))
            sage: Hko = FPModule(Ap, [0], [[Sq(2)], [Sq(1)]])
            sage: f = Hom(Hko, Hko)([(Ap.Sq(0,0,3) + Ap.Sq(0,2,0,1))*Hko.generator(0)])
            sage: f*f == 0
            True
            sage: k = f.kernel_inclusion() # long time
            sage: f.lift(k) # long time
            Module morphism:
              From: Finitely presented left module on 1 generator and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 2, 2, 1]
              To:   Finitely presented left module on 3 generators and 12 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 2, 2, 1]
              Defn: g[0] |--> Sq(1)*g[20]

        Corner cases involving trivial maps::

            sage: M = FPModule(A, [1])
            sage: M1 = FPModule(A, [0])
            sage: M2 = FPModule(A, [0], [[Sq(1)]])
            sage: q = Hom(M1, M2)([M2.generator(0)])
            sage: z = Hom(M, M2).zero()
            sage: lift = z.lift(q)
            sage: lift.domain() is M and lift.codomain() is M1
            True

        .. SEEALSO::

            :meth:`split`
        """
        from sage.modules.free_module_element import vector

        #        self
        #    L -------> N
        #     \         ^
        #       \       |
        #   lift  \     | f
        #           \   |
        #            _| |
        #               M
        L = self.domain()
        N = self.codomain()
        M = f.domain()

        # It is an error to call this function with incompatible arguments.
        if f.codomain() is not N:
            raise ValueError('the codomains of this homomorphism and the homomorphism '
                             'we are lifting over are different')

        # The trivial map lifts over any other map.
        if self.is_zero():
            return Hom(L, M).zero()

        # A non-trivial map never lifts over the trivial map.
        if f.is_zero():
            if verbose:
                print('This homomorphism cannot lift over a trivial homomorphism'
                      ' since it is non-trivial.')
            return None

        xs = [f.solve(self(g)) for g in L.generators()]

        # If some of the generators are not in the image of f, there is no
        # hope finding a lift.
        if None in xs:
            if verbose:
                print('The generators of the domain of this homomorphism do '
                      'not map into the image of the homomorphism we are lifting over.')
            return None

        # If L is free there are no relations to take into consideration.
        if not L.has_relations():
            return Hom(L, M)(xs)

        # The degree of the lifted map f_.
        lift_deg = self.degree() - f.degree()

        # Compute the kernel of f.  The equations we will solve will live in
        # this submodule.
        iK = f.kernel_inclusion(top_dim=max([r.degree() + lift_deg for r in L.relations()]))

        source_degs = [g.degree() + lift_deg for g in L.generators()]
        target_degs = [r.degree() + lift_deg for r in L.relations()]

        # Act on the liftings xs by the relations.
        ys = []
        K = iK.domain()
        all_zero = True

        for r in L.relations():
            target_degree = r.degree() + lift_deg

            y = iK.solve(sum([c * x for c, x in zip(r.dense_coefficient_list(), xs)]))
            if y is None:
                if verbose:
                    print('The homomorphism cannot be lifted in any '
                          'way such that the relations of the domain are '
                          'respected.')
                return None

            if y.is_zero():
                dim = K.vector_presentation(target_degree).dimension()
                ys += dim * [0]  # The zero vector of the appropriate dimension.
            else:
                all_zero = False
                ys += list(y.vector_presentation())

        # If the initial guess already fits the relations, we are done.
        if all_zero:
            return Hom(L, M)(xs)

        block_matrix, R = _create_relations_matrix(
            K, [r.dense_coefficient_list() for r in L.relations()], source_degs, target_degs)

        try:
            solution = R.solve_right(vector(ys))
        except ValueError as error:
            if str(error) == 'matrix equation has no solutions':
                if verbose:
                    print('The homomorphism cannot be lifted in any '
                          'way such that the relations of the domain '
                          'are respected: %s' % error)

                return None
            else:
                raise ValueError(error)

        # Interpret the solution vector as a vector in the direct sum
        # $ K_1\oplus K_2\oplus \ldots \oplus K_n $.
        n = 0
        for j, source_degree in enumerate(source_degs):

            source_dimension = block_matrix[0][j].domain().dimension()

            w = K.element_from_coordinates(
                solution[n:n + source_dimension], source_degree)

            # Subtract the solution w_i from our initial choice of lift
            # for the generator g_i.
            xs[j] -= iK(w)

            n += source_degree

        return Hom(L, M)(xs)

    def split(self, verbose=False):
        r"""
        Return a split of ``self``.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); enable progress messages

        OUTPUT:

        A homomorphism with the property that the composite homomorphism
        `S \circ f = id`, where `S` is ``self``, is The identity homomorphism
        If no such split exist, ``None`` is returned.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FPModule(A, [0,0], [[0, Sq(1)]])
            sage: N = FPModule(A, [0], [[Sq(1)]])
            sage: p = Hom(M, N)([N.generator(0), N.generator(0)])
            sage: s = p.split(); s
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> g[0, 1]
            sage: # Verify that `s` is a splitting:
            sage: p*s
            Module endomorphism of Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> g[0]

        TESTS::

            sage: F = FPModule(A, [0])
            sage: N = FPModule(A, [0], [[Sq(1)]])
            sage: p = Hom(F, N)([N.generator(0)])
            sage: p.split(verbose=True) is None
            The homomorphism cannot be lifted in any way such that the relations
             of the domain are respected: matrix equation has no solutions
            True

        .. SEEALSO::

            :meth:`lift`
        """
        one = End(self.codomain()).identity()
        return one.lift(self, verbose)

    def homology(self, f, top_dim=None, verbose=False):
        r"""
        Compute the sub-quotient module of `H(self, f) =
        \ker(self)/\operatorname{im}(f)` in a range of degrees.

        For a pair of composable morphisms `f: M\to N` and `g: N \to Q` of
        finitely presented modules, the homology module is a finitely
        presented quotient of the kernel sub module `\ker(g) \subset N`.

        INPUT:

        - ``f`` -- a homomorphism with codomain equal to the domain of ``self``
          and image contained in the kernel of this homomorphism
        - ``top_dim`` -- integer (optional); used by this function to stop the
          computation at the given degree
        - ``verbose`` -- boolean (default: ``False``); enable progress messages

        OUTPUT:

        A quotient homomorphism `\ker(self) \to H`, where `H` is isomorphic
        to `H(self, f)` in degrees less than or equal to ``top_dim``.

        .. NOTE::

            If the algebra for this module is finite, then no ``top_dim`` needs
            to be specified in order to ensure that this function terminates.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FPModule(A, [0], [[Sq(3)]])
            sage: N = FPModule(A, [0], [[Sq(2,2)]])
            sage: F = FPModule(A, [0])
            sage: f = Hom(M,N)([A.Sq(2)*N.generator(0)])
            sage: g = Hom(F, M)([A.Sq(4)*A.Sq(1,2)*M.generator(0)])
            sage: ho = f.homology(g)
            sage: ho.codomain()
            Finitely presented left module on 1 generator and 5 relations over
             sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
            sage: ho.codomain().is_trivial()
            False
        """
        k = self.kernel_inclusion(top_dim, verbose)
        f_ = f.lift(k)
        if f_ is None:
            raise ValueError('the image of the given homomorphism is not contained '
                             'in the kernel of this homomorphism; the homology is '
                             'therefore not defined for this pair of maps')

        return f_.cokernel_projection()

    def suspension(self, t):
        r"""
        The suspension of this morphism by the given degree ``t``.

        INPUT:

        - ``t`` -- integer by which the morphism is suspended

        OUTPUT: the morphism which is the suspension of ``self`` by the degree ``t``

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: F1 = FPModule(A, [4,5])
            sage: F2 = FPModule(A, [5])

            sage: f = Hom(F1, F2)( ( F2([Sq(4)]), F2([Sq(5)]) ) ); f
            Module morphism:
              From: Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[4] |--> Sq(4)*g[5]
                    g[5] |--> Sq(5)*g[5]

            sage: e1 = F1([1, 0])
            sage: e2 = F1([0, 1])
            sage: f(e1)
            Sq(4)*g[5]
            sage: f(e2)
            Sq(5)*g[5]

            sage: sf = f.suspension(4); sf
            Module morphism:
              From: Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[8] |--> Sq(4)*g[9]
                    g[9] |--> Sq(5)*g[9]

            sage: sf.domain() is f.domain().suspension(4)
            True

            sage: sf.codomain() is f.codomain().suspension(4)
            True
        """
        if t == 0:
            return self

        D = self.domain().suspension(t)
        C = self.codomain().suspension(t)
        return Hom(D, C)([C(x.lift_to_free().dense_coefficient_list())
                          for x in self._values])

    def cokernel_projection(self):
        r"""
        Return the map to the cokernel of ``self``.

        OUTPUT:

        The natural projection from the codomain of this homomorphism
        to its cokernel.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A1 = SteenrodAlgebra(2, profile=(2,1))
            sage: M = FPModule(A1, [0], [[Sq(2)]])
            sage: F = FPModule(A1, [0])

            sage: r = Hom(F, M)([A1.Sq(1)*M.generator(0)])
            sage: co = r.cokernel_projection(); co
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]
              To:   Finitely presented left module on 1 generator and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]
              Defn: g[0] |--> g[0]

            sage: co.domain().is_trivial()
            False
        """
        new_relations = ([x.dense_coefficient_list()
                          for x in self.codomain().relations()] +
                         [x.dense_coefficient_list() for x in self._values])

        try:
            FPModule = self.base_ring()._fp_graded_module_class
        except AttributeError:
            from .module import FPModule

        coker = FPModule(self.base_ring(),
                         self.codomain().generator_degrees(),
                         relations=tuple(new_relations))

        projection = Hom(self.codomain(), coker)(coker.generators())

        return projection

    def kernel_inclusion(self, top_dim=None, verbose=False):
        r"""
        Return the kernel of ``self``.

        INPUT:

        - ``top_dim`` -- integer (optional); used by this function to stop the
          computation at the given degree
        - ``verbose`` -- boolean (default: ``False``); enable progress messages

        OUTPUT:

        A homomorphism into `\ker(self)` which is an isomorphism in
        degrees less than or equal to ``top_dim``.

        .. NOTE::

            If the algebra for this module is finite, then no ``top_dim`` needs
            to be specified in order to ensure that this function terminates.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: F = FPModule(A3, [1,3]);
            sage: L = FPModule(A3, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]]);
            sage: H = Hom(F, L);

            sage: H([L((A3.Sq(1), 1)), L((0, A3.Sq(2)))]).kernel_inclusion() # long time
            Module morphism:
              From: Finitely presented left module on 2 generators and 1 relation over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [4, 3, 2, 1]
              To:   Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [4, 3, 2, 1]
              Defn: g[3] |--> g[3]
                    g[4] |--> Sq(0,1)*g[1]

            sage: M = FPModule(A3, [0,7], [[Sq(1), 0], [Sq(2), 0], [Sq(4), 0], [Sq(8), Sq(1)], [0, Sq(7)], [0, Sq(0,1,1)+Sq(4,2)]])
            sage: F2 = FPModule(A3, [0], [[Sq(1)], [Sq(2)], [Sq(4)], [Sq(8)], [Sq(15)]])
            sage: H = Hom(M, F2)
            sage: f = H([F2([1]), F2([0])])

            sage: K = f.kernel_inclusion(verbose=True, top_dim=17)
            1. Computing the generators of the kernel presentation:
            Resolving the kernel in the range of dimensions [0, 17]:
             0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17.
            2. Computing the relations of the kernel presentation:
            Computing using the profile:
            (4, 3, 2, 1)
            Resolving the kernel in the range of dimensions [7, 17]:
             7 8 9 10 11 12 13 14 15 16 17.

            sage: K.domain().generators()
            (g[7],)
            sage: K.domain().relations()
            ((Sq(0,1)+Sq(3))*g[7],
             (Sq(0,0,1)+Sq(1,2)+Sq(4,1))*g[7],
             Sq(9)*g[7],
             (Sq(0,1,1)+Sq(4,2))*g[7])

            sage: K
            Module morphism:
              From: Finitely presented left module on 1 generator and 4 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [4, 3, 2, 1]
              To:   Finitely presented left module on 2 generators and 6 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [4, 3, 2, 1]
              Defn: g[7] |--> g[7]
        """
        if verbose:
            print('1. Computing the generators of the kernel presentation:')
        j0 = self._resolve_kernel(top_dim, verbose)
        if verbose:
            print('2. Computing the relations of the kernel presentation:')
        j1 = j0._resolve_kernel(top_dim, verbose)

        # Create a module isomorphic to the ker(self).
        K = j1.fp_module()

        # Return an injection of K into the domain of self such that
        # its image equals ker(self).
        return Hom(K, j0.codomain())(j0._values)

    def image(self, top_dim=None, verbose=False):
        r"""
        Compute the image of ``self``.

        INPUT:

        - ``top_dim`` -- integer (optional); used by this function to stop the
          computation at the given degree
        - ``verbose`` -- boolean (default: ``False``); enable progress messages

        OUTPUT:

        A homomorphism into `\operatorname{im}(self)` that is an
        isomorphism in degrees less than or equal to ``top_dim``

        .. NOTE::

            If the algebra for this module is finite, then no ``top_dim``
            needs to be specified in order to ensure that this function
            terminates.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A3 = SteenrodAlgebra(2, profile=(4,3,2,1))
            sage: F = FPModule(A3, [1,3]);
            sage: L = FPModule(A3, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]]);
            sage: H = Hom(F, L);

            sage: H([L((A3.Sq(1), 1)), L((0, A3.Sq(2)))]).image() # long time
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [4, 3, 2, 1]
              To:   Finitely presented left module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [4, 3, 2, 1]
              Defn: g[3] |--> Sq(1)*g[2] + g[3]

            sage: M = FPModule(A3, [0,7], [[Sq(1), 0], [Sq(2), 0], [Sq(4), 0], [Sq(8), Sq(1)], [0, Sq(7)], [0, Sq(0,1,1)+Sq(4,2)]])
            sage: F2 = FPModule(A3, [0], [[Sq(1)], [Sq(2)], [Sq(4)], [Sq(8)], [Sq(15)]])
            sage: H = Hom(M, F2)
            sage: f = H([F2([1]), F2([0])])
            sage: K = f.image(verbose=True, top_dim=17)
            1. Computing the generators of the image presentation:
            Resolving the image in the range of dimensions [0, 17]:
             0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17.
            2. Computing the relations of the image presentation:
            Computing using the profile:
            (4, 3, 2, 1)
            Resolving the kernel in the range of dimensions [0, 17]:
             0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17.

            sage: K.is_injective()  # long time
            True
            sage: K.domain().generator_degrees()
            (0,)
            sage: K.domain().relations()
            (Sq(1)*g[0], Sq(2)*g[0], Sq(4)*g[0], Sq(8)*g[0])
            sage: K.domain().is_trivial()
            False
        """
        if verbose:
            print('1. Computing the generators of the image presentation:')
        j0 = self._resolve_image(top_dim, verbose)
        if verbose:
            print('2. Computing the relations of the image presentation:')
        j1 = j0._resolve_kernel(top_dim, verbose)

        # Create a module isomorphic to the im(self).
        I = j1.fp_module()

        # Return an injection of I into the codomain of self such that
        # its image equals im(self)
        return Hom(I, j0.codomain())(j0._values)

    def is_injective(self, top_dim=None, verbose=False):
        r"""
        Return ``True`` if and only if ``self`` has a trivial kernel.

        INPUT:

        - ``top_dim`` -- integer (optional); used by this function to stop the
          computation at the given degree
        - ``verbose`` -- boolean (default: ``False``); enable progress messages

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)

            sage: K = FPModule(A, [2], [[Sq(2)]])
            sage: HZ = FPModule(A, [0], [[Sq(1)]])

            sage: f = Hom(K, HZ)([Sq(2)*HZ([1])])
            sage: f.is_injective(top_dim=23)
            True

        TESTS::

            sage: Z = FPModule(A, [])
            sage: Hom(Z, HZ).zero().is_injective(top_dim=8)
            True
        """
        j0 = self._resolve_kernel(top_dim, verbose)
        return j0.domain().is_trivial()

    def is_surjective(self):
        r"""
        Return ``True`` if and only if ``self`` has a trivial cokernel.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: F = FPModule(A, [0])

            sage: f = Hom(F,F)([Sq(1)*F.generator(0)])
            sage: f.is_surjective()
            False

        TESTS::

            sage: Z = FPModule(A, [])
            sage: Hom(F, Z).zero().is_surjective()
            True
        """
        return self.cokernel_projection().is_zero()

    def _resolve_kernel(self, top_dim=None, verbose=False):
        r"""
        Resolve the kernel of ``self`` by a free module.

        INPUT:

        - ``top_dim`` -- integer (optional); used by this function to stop the
          computation at the given degree
        - ``verbose`` -- boolean (default: ``False``); enable progress messages

        OUTPUT:

        A homomorphism `j: F \rightarrow D`, where `D` is the domain of
        ``self`` and `F` is free and such that `\ker(self) = \operatorname{im}(j)`
        in all degrees less than or equal to ``top_dim``.

        .. NOTE::

            If the algebra for this module is finite dimensional, then no
            ``top_dim`` needs to be specified in order to ensure that
            this function terminates.

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: s = SymmetricFunctions(QQ).s()
            sage: F = s.free_graded_module([0,0])
            sage: L = FPModule(s, [0,0], [[s[3],s[2,1]], [0,s[2]]])
            sage: f = Hom(F, L)([L([s[2], 0]), L([0, s[2]])])
            sage: f._resolve_kernel()  # long time
            Traceback (most recent call last):
            ...
            ValueError: a top dimension must be specified for this calculation to terminate
            sage: f._resolve_kernel(top_dim=10)
            Module morphism:
              From: Free graded left module on 2 generators over Symmetric Functions over Rational Field in the Schur basis
              To:   Free graded left module on 2 generators over Symmetric Functions over Rational Field in the Schur basis
              Defn: s[]*g[0] |--> s[]*g[0, 1]
                    s[]*g[3] |--> s[3]*g[0, 0]
        """
        # Let
        #
        #  1) `j` be a homomorphism into `\ker(self)`, and
        #  2) 'n' be an integer.
        #
        # The induction loop starts each iteration assuming that `j` is onto
        # the kernel in degrees below `n`.  Each iteration of the loop then
        # extends the map `j` minimally so that `j_n` becomes onto the kernel.
        #
        # This induction step is then repeated for all `n \leq` ``top_dim``.
        #

        domain = self.domain()
        R = self.base_ring()

        if self.is_zero():
            # Epsilon: F_0 -> M
            F_0 = R.free_graded_module(domain.generator_degrees())
            epsilon = Hom(F_0, domain)(tuple(domain.generators()))
            return epsilon

        # Create the trivial module F_ to start with.
        F_ = R.free_graded_module(())
        j = Hom(F_, domain).zero()

        dim = domain.connectivity()
        if dim == infinity:
            if verbose:
                print('The domain of the morphism is trivial, so there is nothing to resolve.')
            return j

        if not R.dimension() < infinity:
            limit = infinity
        else:
            limit = _top_dim(R) + max(domain.generator_degrees())

        if top_dim is not None:
            limit = min(top_dim, limit)

        if limit == infinity:
            raise ValueError('a top dimension must be specified for this calculation to terminate')

        if verbose:
            if dim > limit:
                print('The dimension range is empty: [%d, %d]' % (dim, limit))
            else:
                print('Resolving the kernel in the range of dimensions [%d, %d]:' % (dim, limit), end='')

        # The induction loop.
        for n in range(dim, limit + 1):

            if verbose:
                print(' %d' % n, end='')

            # We have taken care of the case when self is zero, so the
            # vector presentation exists.
            self_n = self.vector_presentation(n)
            kernel_n = self_n.kernel()

            if not kernel_n:
                continue

            generator_degrees = tuple(x.degree() for x in F_.generators())

            if j.is_zero():
                # The map j is not onto in degree `n` of the kernel.
                new_generator_degrees = kernel_n.rank() * (n,)
                F_ = R.free_graded_module(generator_degrees + new_generator_degrees)

                new_values = tuple([
                    domain.element_from_coordinates(q, n) for q in kernel_n.basis()])

            else:
                Q_n = kernel_n.quotient(j.vector_presentation(n).image())

                if not Q_n.rank():
                    continue

                # The map j is not onto in degree `n` of the kernel.
                new_generator_degrees = Q_n.rank() * (n,)
                F_ = R.free_graded_module(generator_degrees + new_generator_degrees)

                new_values = tuple([
                    domain.element_from_coordinates(Q_n.lift(q), n) for q in Q_n.basis()])

            # Create a new homomorphism which is surjective onto the kernel
            # in all degrees less than, and including `n`.
            j = Hom(F_, domain)(j._values + new_values)

        if verbose:
            print('.')
        return j

    def _resolve_image(self, top_dim=None, verbose=False):
        r"""
        Resolve the image of ``self`` by a free module.

        INPUT:

        - ``top_dim`` -- integer (optional); used by this function to stop the
          computation at the given degree
        - ``verbose`` -- boolean (default: ``False``); enable progress messages

        OUTPUT:

        A homomorphism `j: F \rightarrow C` where `C` is the codomain
        of ``self``, where `F` is free, and
        `\operatorname{im}(self) = \operatorname{im}(j)` in all degrees less
        than or equal to ``top_dim``.

        .. NOTE::

            If the algebra for this module is finite, then no ``top_dim`` needs
            to be specified in order to ensure that this function terminates.

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: s = SymmetricFunctions(QQ).s()
            sage: F = s.free_graded_module([0,0])
            sage: L = FPModule(s, [0,0], [[s[3],s[2,1]], [0,s[2]]])
            sage: f = Hom(F, L)([L([s[2], 0]), L([0, s[2]])])
            sage: f._resolve_image()
            Traceback (most recent call last):
            ...
            ValueError: a top dimension must be specified for this calculation to terminate
            sage: f._resolve_image(top_dim=10)
            Module morphism:
              From: Free graded left module on 1 generator over Symmetric Functions over Rational Field in the Schur basis
              To:   Finitely presented left module on 2 generators and 2 relations over Symmetric Functions over Rational Field in the Schur basis
              Defn: s[]*g[2] |--> s[2]*g[0, 0]
        """
        # Let
        #
        #  1) `j` be a homomorphism into `\im(self)`, and
        #  2) 'n' be an integer.
        #
        # The induction loop starts each iteration assuming that `j` is onto
        # the image in degrees below `n`.  Each iteration of the loop then
        # extends the map `j` minimally so that `j_n` becomes onto the image.
        #
        # This induction step is then repeated for all `n \leq` ``top_dim``.
        #

        # Create the trivial module F_ to start with.
        R = self.base_ring()
        codomain = self.codomain()
        F_ = R.free_graded_module(())
        j = Hom(F_, codomain).zero()

        dim = self.codomain().connectivity()
        if dim == infinity:
            if verbose:
                print('The codomain of the morphism is trivial, so there is nothing to resolve.')
            return j

        try:
            self_degree = self.degree()
        except ValueError:
            if verbose:
                print('The homomorphism is trivial, so there is nothing to resolve.')
            return j

        degree_values = [0] + [v.degree() for v in self._values if v]
        limit = (infinity if not R.dimension() < infinity else
                 (_top_dim(R) + max(degree_values)))

        if top_dim is not None:
            limit = min(top_dim, limit)

        if limit == infinity:
            raise ValueError('a top dimension must be specified for this calculation to terminate')

        if verbose:
            if dim > limit:
                print('The dimension range is empty: [%d, %d]' % (dim, limit))
            else:
                print('Resolving the image in the range of dimensions [%d, %d]:' % (dim, limit), end='')

        for n in range(dim, limit + 1):

            if verbose:
                print(' %d' % n, end='')

            self_n = self.vector_presentation(n - self_degree)
            image_n = self_n.image()

            if image_n.dimension() == 0:
                continue

            generator_degrees = tuple(x.degree() for x in F_.generators())
            if j.is_zero():
                # The map j is not onto in degree `n` of the image.
                new_generator_degrees = image_n.rank() * (n,)
                F_ = R.free_graded_module(generator_degrees + new_generator_degrees)

                new_values = tuple([
                    self.codomain().element_from_coordinates(q, n) for q in image_n.basis()])

            else:

                j_n = j.vector_presentation(n)
                Q_n = image_n.quotient(j_n.image())

                if Q_n.dimension() == 0:
                    continue

                # The map j is not onto in degree `n` of the image.
                new_generator_degrees = Q_n.rank() * (n,)
                F_ = R.free_graded_module(generator_degrees + new_generator_degrees)

                new_values = tuple([
                    self.codomain().element_from_coordinates(Q_n.lift(q), n) for q in Q_n.basis()])

            # Create a new homomorphism which is surjective onto the image
            # in all degrees less than, and including `n`.
            j = Hom(F_, self.codomain())(j._values + new_values)

        if verbose:
            print('.')
        return j

    def fp_module(self):
        r"""
        Create a finitely presented module from ``self``.

        OUTPUT:

        The finitely presented module having presentation equal to ``self``
        as long as the domain and codomain are free.

        EXAMPLES:

        We construct examples with free modules that are presented
        with a redundant relation::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: F1 = FPModule(A, (2,), [[0]])
            sage: F2 = FPModule(A, (0,), [[0]])
            sage: v = F2([Sq(2)])
            sage: pres = Hom(F1, F2)([v])
            sage: M = pres.fp_module(); M
            Finitely presented left module on 1 generator and 1 relation over
             mod 2 Steenrod algebra, milnor basis
            sage: M.generator_degrees()
            (0,)
            sage: M.relations()
            (Sq(2)*g[0],)

            sage: F2 = A.free_graded_module((0,))
            sage: v = F2([Sq(2)])
            sage: pres = Hom(F1, F2)([v])
            sage: M = pres.fp_module(); M
            Finitely presented left module on 1 generator and 1 relation over
             mod 2 Steenrod algebra, milnor basis
            sage: M.generator_degrees()
            (0,)
            sage: M.relations()
            (Sq(2)*g[0],)

            sage: F3 = FPModule(A, (0,), [[Sq(4)]])
            sage: v = F3([Sq(2)])
            sage: pres = Hom(F1, F3)([v])
            sage: pres.fp_module()
            Traceback (most recent call last):
            ...
            ValueError: this is not a morphism between free modules
        """
        if self.domain().has_relations() or self.codomain().has_relations():
            raise ValueError("this is not a morphism between free modules")
        try:
            FPModule = self.base_ring()._fp_graded_module_class
        except AttributeError:
            from .module import FPModule
        return FPModule(self.base_ring(),
                        self.codomain().generator_degrees(),
                        tuple([r.dense_coefficient_list() for r in self._values]))


@cached_function
def _top_dim(algebra):
    r"""
    The top dimension of ``algebra``.

    This returns infinity if the algebra is infinite-dimensional. If
    the algebra has a ``top_class`` method, then it is used in the
    computation; otherwise the computation may be slow.

    EXAMPLES::

        sage: from sage.modules.fp_graded.morphism import _top_dim
        sage: E.<x,y,z> = ExteriorAlgebra(QQ)
        sage: _top_dim(E)
        3

        sage: _top_dim(SteenrodAlgebra(2, profile=(3,2,1)))
        23
    """
    if not algebra.dimension() < infinity:
        return infinity
    try:
        alg_top_dim = algebra.top_class().degree()
    except AttributeError:
        # This could be very slow, but it should handle the case
        # when the algebra is finite-dimensional but has no
        # top_class method, e.g., exterior algebras.
        alg_top_dim = max(a.degree() for a in algebra.basis())
    return alg_top_dim
