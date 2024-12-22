r"""
Finitely presented graded modules over the Steenrod algebra

This package allows the user to define finitely presented modules
over the mod `p` Steenrod algebra, elements of them, and
morphisms between them.  Methods are provided for doing basic homological
algebra, e.g. computing kernels and images of homomorphisms, and finding
free resolutions of modules.

.. RUBRIC:: Theoretical background

The category of finitely presented graded modules over an arbitrary
non-Noetherian graded ring `R` is not abelian in general, since
kernels of homomorphisms are not necessarily finitely presented.

The mod `p` Steenrod algebra is non-Noetherian, but it is the
union of a countable set of finite sub-Hopf algebras
([Mar1983]_ Ch. 15, Sect. 1, Prop 7). It is therefore an example of a
`P`-algebra ([Mar1983]_ Ch. 13).

Any finitely presented module over the Steenrod algebra is defined
over one of these finite sub-Hopf algebras.  Similarly, any
homomorphism between finitely presented modules over the Steenrod
algebra is defined over a finite sub-Hopf algebra of the Steenrod
algebra.  Further, tensoring up from the category of modules over a
sub-Hopf algebra is an exact functor, since the Steenrod algebra is
free over any sub-Hopf algebra.

It follows that kernels, cokernels, images, and, more generally, any finite
limits or colimits can be computed in the category of modules over the
Steenrod algebra by working in the category of modules over an appropriate
finite sub-Hopf algebra.

It is also the case that presentations of modules and the images of the
generators of the domain of a homomorphism are the same over the sub-Hopf
algebra and over the whole Steenrod algebra, so that the tensoring up is
entirely implicit and requires no computation.

The definitions and computations carried out by this package are thus done
as follows.   First, find a small finite sub-Hopf algebra over which the
computation can be done.   Then, carry out the calculation there, where it
is a finite problem and can be reduced to linear algebra over a finite
prime field.

For examples, see the `Steenrod algebra modules
<../../../../../../thematic_tutorials/steenrod_algebra_modules.html>`_
thematic tutorial.

TESTS::

    sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
    sage: A = SteenrodAlgebra(2)
    sage: K = SteenrodFPModule(A, [1,3]); K
    Free graded left module on 2 generators over ...
    sage: K.category()
    Category of finite dimensional graded modules with basis over mod 2 Steenrod algebra, milnor basis
    sage: L = SteenrodFPModule(A, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]]); L
    Finitely presented left module on 2 generators and 2 relations ...
    sage: M = SteenrodFPModule(A, [2,3], [[Sq(2),Sq(1)]]); M
    Finitely presented left module on 2 generators and 1 relation ...
    sage: m = M((0,1)); m
    g[3]
    sage: K.is_parent_of(m)
    False
    sage: L.is_parent_of(m)
    False
    sage: M.is_parent_of(m)
    True

    sage: SteenrodFPModule(ZZ, [0])
    Traceback (most recent call last):
    ...
    AttributeError: 'sage.rings.integer_ring.IntegerRing_class' object has no attribute 'free_graded_module'...

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

# ****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu>
#             and          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.infinity import infinity
from sage.algebras.steenrod.steenrod_algebra import SteenrodAlgebra_generic
from sage.modules.fp_graded.module import FPModule
from sage.modules.fp_graded.free_module import FreeGradedModule
from .profile import enveloping_profile_elements


class SteenrodModuleMixin:
    """
    Mixin class for common methods of the Steenrod algebra modules.
    """
    def profile(self):
        r"""
        Return a finite profile over which ``self`` can be defined.

        Any finitely presented module over the Steenrod algebra can be
        defined over a finite-dimensional sub-Hopf algebra, and this
        method identifies such a sub-Hopf algebra and returns its
        profile function.

        .. NOTE::

            The profile produced by this function is reasonably small
            but is not guaranteed to be minimal.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = SteenrodFPModule(A, [0,1], [[Sq(2),Sq(1)],[0,Sq(2)],[Sq(3),0]])
            sage: M.profile()
            (2, 1)

        TESTS::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: X = SteenrodFPModule(A, [0])
            sage: X.profile()
            (0,)
        """
        elements = [coeffifient for value in self.relations()
                    for coeffifient in value.dense_coefficient_list()]
        elements = [a for a in elements if a not in (0, 1)]

        return enveloping_profile_elements(elements,
                                           char=self.base_ring().characteristic())

    def export_module_definition(self, powers_of_two_only=True):
        r"""
        Return the module to the input
        `format used by R. Bruner's Ext software
        <http://www.math.wayne.edu/~rrb/cohom/modfmt.html>`_ as a string.

        INPUT:

        - ``powers_of_two_only`` -- boolean (default: ``True``); if the
          output should contain the action of all Steenrod squaring operations
          (restricted by the profile), or only the action of the operations
          of degree equal to a power of two

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A1 = algebra=SteenrodAlgebra(p=2, profile=[2,1])
            sage: M = SteenrodFPModule(A1, [0])
            sage: print(M.export_module_definition())
            8 0 1 2 3 3 4 5 6
            0 1 1 1
            2 1 1 4
            3 1 1 5
            6 1 1 7
            0 2 1 2
            1 2 2 3 4
            2 2 1 5
            3 2 1 6
            4 2 1 6
            5 2 1 7
            sage: N = SteenrodFPModule(A1, [0], [[Sq(1)]])
            sage: print(N.export_module_definition())
            4 0 2 3 5
            1 1 1 2
            0 2 1 1
            2 2 1 3
            sage: print(N.export_module_definition(powers_of_two_only=False))
            4 0 2 3 5
            1 1 1 2
            0 2 1 1
            2 2 1 3
            0 3 1 2
            sage: A2 = SteenrodAlgebra(p=2, profile=[3,2,1])
            sage: Hko = SteenrodFPModule(A2, [0], [[Sq(1)], [Sq(2)]])
            sage: print(Hko.export_module_definition())
            8 0 4 6 7 10 11 13 17
            2 1 1 3
            4 1 1 5
            1 2 1 2
            5 2 1 6
            0 4 1 1
            2 4 1 4
            3 4 1 5
            6 4 1 7

        TESTS::

            sage: A = SteenrodAlgebra()
            sage: M = A.free_graded_module([])
            sage: M.export_module_definition()
            Traceback (most recent call last):
            ...
            ValueError: this module is not defined over a finite algebra
            sage: A1 = SteenrodAlgebra(profile=[2,1])
            sage: M1 = A1.free_graded_module([])
            sage: s = M1.export_module_definition()
            The module connectivity is infinite, so there is nothing to export.
            sage: s
            ''
            sage: P1 = SteenrodAlgebra(p=5, profile=[[], [2,1]])
            sage: N = P1.free_graded_module([])
            sage: N.export_module_definition()
            Traceback (most recent call last):
            ...
            ValueError: this function is not implemented for odd primes
        """
        if not self.base_ring().is_finite():
            raise ValueError('this module is not defined over a finite algebra')
            return

        if self.base_ring().characteristic() != 2:
            raise ValueError('this function is not implemented for odd primes')
            return

        n = self.connectivity()
        if n == infinity:
            print('The module connectivity is infinite, so there is ' +
                  'nothing to export.')
            return ''

        limit = self.base_ring().top_class().degree() + max(self.generator_degrees())

        # Create a list of bases, one for every module degree we consider.
        vector_space_basis = [self.basis_elements(i)
                              for i in range(n, limit + 1)]

        additive_generator_degrees = []
        additive_generator_global_indices = [0]
        for dim, basis_vectors in enumerate(vector_space_basis):
            additive_generator_global_indices.append(
                len(basis_vectors) + additive_generator_global_indices[-1])
            additive_generator_degrees += len(basis_vectors) * [dim + n]

        # Print the degrees of the additive generators.
        ret = '%d %s' % (len(additive_generator_degrees),
                         ' '.join('%d' % x for x in additive_generator_degrees))

        # A private function which transforms a vector in a given dimension
        # to a vector of global indices for the basis elements corresponding
        # to the nonzero entries in the vector.  E.g.
        # _GetIndices(dim=2, vec=(1,0,1)) will return a vector of length two,
        # (a, b), where a is the index of the first vector in the basis for
        # the 2-dimensional part of the module, and b is the index of the
        # last vector in the same part.
        def _GetIndices(dim, vec):
            if len(vector_space_basis[dim]) != len(vec):
                raise ValueError('the given vector\n%s\nhas the wrong size, it should be %d' % (str(vec), len(vector_space_basis[dim])))
            base_index = additive_generator_global_indices[dim]
            return [base_index + a for a, c in enumerate(vec) if c != 0]

        profile = self.base_ring()._profile
        if powers_of_two_only:
            powers = [2**i for i in range(profile[0])]
        else:
            powers = range(1, 2**profile[0])

        R = self.base_ring()
        for k in powers:
            Sqk = R.Sq(k)
            images = [[(Sqk * x).vector_presentation() for x in D]
                      for D in vector_space_basis]

            element_index = 0

            # Note that the variable dim is relative to the bottom dimension, n.
            for dim, image in enumerate(images):
                for im in image:
                    if im != 0 and im is not None:
                        values = _GetIndices(dim + k, im)

                        ret += "\n%d %d %d %s" % (
                            element_index,
                            k,
                            len(values),
                            " ".join("%d" % x for x in values))
                    element_index += 1
        return ret


class SteenrodFPModule(FPModule, SteenrodModuleMixin):
    r"""
    Create a finitely presented module over the Steenrod algebra.

    .. SEEALSO::

        The thematic tutorial on `Steenrod algebra modules
        <../../../../../../thematic_tutorials/steenrod_algebra_modules.html>`_.

    INPUT:

    One of the following:

    - ``arg0`` -- a morphism such that the module is the cokernel, or
      a free graded module, in which case the output is the same
      module, viewed as finitely presented

    Otherwise:

    - ``arg0`` -- the graded connected algebra over which the module is
      defined; this algebra must be equipped with a graded basis

    - ``generator_degrees`` -- tuple of integer degrees

    - ``relations`` -- tuple of relations; a relation is a tuple of
      coefficients `(c_1, \ldots, c_n)`, ordered so that they
      correspond to the module generators, that is, such a tuple
      corresponds to the relation

      .. MATH::

          c_1 g_1 + \ldots + c_n g_n = 0

      if the generators are `(g_1, \ldots, g_n)`

    TESTS::

        sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
        sage: SteenrodFPModule(SteenrodAlgebra(2), (0,))
        Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
    """
    def _Hom_(self, other, category=None):
        """
        The homset from ``self`` to ``other``.

        TESTS::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: F = SteenrodFPModule(A, [1,3])
            sage: L = SteenrodFPModule(A, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])

            sage: Hom(F, L)
            Set of Morphisms from Free graded left module on 2 generators ...

            sage: M = A.free_graded_module((2,1))
            sage: Hom(F, M)
            Set of Morphisms from Free graded left module on 2 generators ...
        """
        from .homspace import SteenrodFPModuleHomspace
        return SteenrodFPModuleHomspace(self, other, category=category)

    def resolution(self, k, top_dim=None, verbose=False):
        r"""
        A free resolution of ``self`` of the given length.

        INPUT:

        - ``k`` -- nonnegative integer
        - ``top_dim`` -- (optional) stop the computation at this degree
        - ``verbose`` -- boolean (default: ``False``); whether log messages are
          printed

        OUTPUT: list of homomorphisms `[\epsilon, f_1, \ldots, f_k]` such that

        .. MATH::

            \begin{gathered}
            f_i: F_i \to F_{i-1} \text{ for } 1\leq i\leq k, \\
            \epsilon: F_0\to M,
            \end{gathered}

        where each `F_i` is a finitely generated free module, and the
        sequence

        .. MATH::

            F_k \xrightarrow{f_k} F_{k-1} \xrightarrow{f_{k-1}} \ldots
            \rightarrow F_0 \xrightarrow{\epsilon} M \rightarrow 0

        is exact. Note that the 0th element in this list is the map
        `\epsilon`, and the rest of the maps are between free
        modules.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: Hko = SteenrodFPModule(A, [0], [[Sq(1)], [Sq(2)]])

            sage: res = Hko.resolution(5, verbose=True)
            Computing f_1 (1/5)
            Computing f_2 (2/5)
            Computing using the profile:
            (2, 1)
            Resolving the kernel in the range of dimensions [1, 8]: 1 2 3 4 5 6 7 8.
            Computing f_3 (3/5)
            Computing using the profile:
            (2, 1)
            Resolving the kernel in the range of dimensions [2, 10]: 2 3 4 5 6 7 8 9 10.
            Computing f_4 (4/5)
            Computing using the profile:
            (2, 1)
            Resolving the kernel in the range of dimensions [3, 13]: 3 4 5 6 7 8 9 10 11 12 13.
            Computing f_5 (5/5)
            Computing using the profile:
            (2, 1)
            Resolving the kernel in the range of dimensions [4, 18]: 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18.
            sage: [x.domain() for x in res]
            [Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 3 generators over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 4 generators over mod 2 Steenrod algebra, milnor basis]

        When there are no relations, the resolution is trivial::

            sage: M = SteenrodFPModule(A, [0])
            sage: M.resolution(4)
            [Module endomorphism of Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
               Defn: g[0] |--> g[0],
             Module morphism:
               From: Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
               To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis,
             Module endomorphism of Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis,
             Module endomorphism of Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis,
             Module endomorphism of Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis]
        """
        algebra = self.base_ring()
        finite_algebra = SteenrodAlgebra_generic(algebra.prime(),
                                                 profile=self.profile())

        # Change rings to the finite algebra, and call the base class
        # implementation of this function.
        res = FPModule.resolution(self.change_ring(finite_algebra),
                                  k,
                                  top_dim=top_dim,
                                  verbose=verbose)

        # Change rings back to the original Steenrod algebra.
        # Also convert the maps and modules from FPModule to SteenrodFPModule.
        return [j.change_ring(self.base_ring()) for j in res]


class SteenrodFreeModule(FreeGradedModule, SteenrodModuleMixin):
    def _Hom_(self, Y, category):
        r"""
        The internal hook used by the free function
        :meth:`sage.categories.homset.hom.Hom` to create homsets
        involving ``self``.

        TESTS::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: M._Hom_(M, category=None)
            Set of Morphisms from Free graded left module on 2 generators
              over mod 2 Steenrod algebra, milnor basis
             to Free graded left module on 2 generators
              over mod 2 Steenrod algebra, milnor basis
             in Category of finite dimensional graded modules with basis
              over mod 2 Steenrod algebra, milnor basis

            sage: from sage.modules.fp_graded.module import FPModule
            sage: F = FPModule(A, [1,3])
            sage: Hom(M, F)
            Set of Morphisms from Free graded left module on 2 generators ...
        """
        from sage.modules.fp_graded.steenrod.homspace import SteenrodFreeModuleHomspace
        return SteenrodFreeModuleHomspace(self, Y, category)
