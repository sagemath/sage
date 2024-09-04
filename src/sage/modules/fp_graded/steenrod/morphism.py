r"""
Homomorphisms of finitely presented modules over the Steenrod algebra

This class implements construction and basic manipulation of homomorphisms
between :mod:`finitely presented graded modules
<sage.modules.fp_graded.steenrod.module>` over the mod `p` Steenrod algebra.

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

from sage.categories.homset import Hom

from sage.algebras.steenrod.steenrod_algebra import SteenrodAlgebra_generic
from sage.modules.fp_graded.morphism import FPModuleMorphism
from sage.modules.fp_graded.free_morphism import FreeGradedModuleMorphism
from .profile import enveloping_profile_elements


class SteenrodFPModuleMorphism(FPModuleMorphism):
    def profile(self):
        r"""
        Return a finite profile over which ``self`` can be defined.

        This is in some ways the key method for these morphisms. As
        discussed in the "Theoretical background" section of
        :mod:`sage.modules.fp_graded.steenrod.module`, any
        homomorphism of finitely presented modules over the Steenrod
        algebra can be defined over a finite-dimensional sub-Hopf
        algebra, and this method identifies such a sub-Hopf algebra
        and returns its profile function.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = SteenrodFPModule(A, [0,1], [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: one = Hom(M,M).identity()
            sage: one.profile()
            (2, 1)
            sage: zero = Hom(M,M).zero()
            sage: zero.profile()
            (2, 1)
            sage: A_fin = SteenrodAlgebra(2, profile=(2,1))
            sage: M_fin = M.change_ring(A_fin)

        Change the ring of the module ``M``::

            sage: M_fin.change_ring(A) is M
            True

        We can change rings to the finite sub-Hopf algebra defined by
        the profile we just computed::

            sage: one_fin = one.change_ring(A_fin)
            sage: one_fin.domain()
            Finitely presented left module on 2 generators and 2 relations over
             sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]

        If we change back to the full Steenrod algebra, we are back where
        we started::

            sage: one_fin.change_ring(A) == one
            True
        """
        def _flatten(f):
            return [c for value in f for c in value.dense_coefficient_list()]

        elements = (_flatten(self.domain().relations())
                    + _flatten(self.codomain().relations())
                    + _flatten(self.values()))
        elements = [a for a in elements if a not in (0, 1)]

        return enveloping_profile_elements(elements,
                                           char=self.base_ring().characteristic())

    def is_injective(self, top_dim=None, verbose=False):
        r"""
        Return ``True`` if ``self`` is injective.

        INPUT:

        - ``top_dim`` -- (optional) stop the computation at this degree; if
          not specified, this is determined using :meth:`profile`
        - ``verbose`` -- boolean (default: ``False``); whether log messages are
          printed

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = SteenrodFPModule(A, [0,1], [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = SteenrodFPModule(A, [0], [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: g = Hom(S, M).zero()
            sage: g.is_injective()
            False
            sage: z = Hom(SteenrodFPModule(A, []), M).zero()
            sage: z.is_injective()
            True
            sage: z.is_zero()
            True
        """
        algebra = self.base_ring()
        finite_algebra = SteenrodAlgebra_generic(algebra.prime(), profile=self.profile())
        return FPModuleMorphism.is_injective(self.change_ring(finite_algebra),
                                             top_dim=top_dim, verbose=verbose)

    def kernel_inclusion(self, top_dim=None, verbose=False):
        r"""
        Return the kernel of ``self`` as a morphism.

        INPUT:

        - ``top_dim`` -- (optional) stop the computation at this degree; if
          not specified, this is determined using :meth:`profile`
        - ``verbose`` -- boolean (default: ``False``); whether log messages are
          printed

        OUTPUT: an injective homomorphism into the domain ``self`` which is
        onto the kernel of this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = SteenrodFPModule(A, [0,1], [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = SteenrodFPModule(A, [0], [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: k = f.kernel_inclusion()
            sage: k == 0
            True

        Since k is both trivial and injective, its domain should
        be the zero module::

            sage: k.domain().is_trivial()
            True

            sage: g = Hom(S, M)([M([Sq(3),Sq(2)])])
            sage: h = g.kernel_inclusion()
            sage: h.is_identity()
            True
            sage: ker = h.domain();
            sage: ker is S
            True

        So `g` had to be trivial::

            sage: g.is_zero()
            True
        """
        return self._action(FPModuleMorphism.kernel_inclusion, top_dim=top_dim, verbose=verbose)

    def cokernel_projection(self, verbose=False):
        r"""
        Compute the map to the cokernel of ``self``.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); whether log messages are
          printed

        OUTPUT:

        The natural projection from the codomain of this homomorphism
        to its cokernel.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A1 = SteenrodAlgebra(2, profile=(2,1))
            sage: M = SteenrodFPModule(A1, [0], [[Sq(2)]])
            sage: F = SteenrodFPModule(A1, [0])

            sage: r = Hom(F, M)([A1.Sq(1)*M.generator(0)])
            sage: co = r.cokernel_projection(); co
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]
              To:   Finitely presented left module on 1 generator and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]
              Defn: g[0] |--> g[0]

            sage: co.domain().is_trivial()
            False
        """
        from .module import SteenrodFPModule
        new_relations = ([x.dense_coefficient_list()
                          for x in self.codomain().relations()] +
                         [x.dense_coefficient_list() for x in self._values])

        coker = SteenrodFPModule(self.base_ring(),
                                 self.codomain().generator_degrees(),
                                 relations=tuple(new_relations))

        projection = Hom(self.codomain(), coker)(coker.generators())

        return projection

    def image(self, top_dim=None, verbose=False):
        r"""
        Return the image of ``self``.

        INPUT:

        - ``top_dim`` -- integer (optional); used by this function to stop the
          computation at the given degree
        - ``verbose`` -- boolean (default: ``False``); whether log messages are
          printed

        OUTPUT:

        An injective homomorphism into the codomain of ``self`` which is
        onto the image of ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: M = SteenrodFPModule(A, [0,1], [[Sq(2),Sq(1)], [0,Sq(2)]])
            sage: S = SteenrodFPModule(A, [0], [[Sq(2)]])
            sage: f = Hom(S, M)([M([0,1])])
            sage: f.is_injective()
            True
            sage: i = f.image(); i
            Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[1] |--> g[1]
            sage: i.codomain() is M
            True

        Lift the map ``f`` over the inclusion ``i``::

            sage: f_ = f.lift(i)
            sage: f_.is_injective()
            True
            sage: f_.is_surjective()
            True

            sage: g = Hom(S, M)([M([Sq(3),Sq(2)])])
            sage: j = g.image(); j
            Module morphism:
              From: Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 2 relations over mod 2 Steenrod algebra, milnor basis

        So ``g`` had to be trivial::

            sage: g.is_zero()
            True
        """
        return self._action(FPModuleMorphism.image, top_dim=top_dim, verbose=verbose)

    def _resolve_kernel(self, top_dim=None, verbose=False):
        r"""
        Resolve the kernel of this homomorphism by a free module.

        INPUT:

        - ``top_dim`` -- (optional) stop the computation at this degree; if
          not specified, this is determined using :meth:`profile`
        - ``verbose`` -- boolean (default: ``False``); whether log messages are
          printed

        OUTPUT: a homomorphism `j: F \rightarrow D` where `D` is the domain of
        this homomorphism, `F` is free and such that `\ker(self) = \operatorname{im}(j)`.

        TESTS::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: F = SteenrodFPModule(A, [0,0])
            sage: L = SteenrodFPModule(A, [0,0], [[Sq(3),Sq(0,1)], [0,Sq(2)]])
            sage: f = Hom(F, L)([L([Sq(2),0]), L([0, Sq(2)])])
            sage: f._resolve_kernel()
            Module morphism:
              From: Free graded left module on 3 generators over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
              Defn: g[0, 0] |--> g[0, 1]
                    g[3, 0] |--> Sq(0,1)*g[0, 0]
                    g[3, 1] |--> Sq(3)*g[0, 0]

        An odd primary example::

            sage: A3 = SteenrodAlgebra(3)
            sage: F0 = A3.free_graded_module([32, 40])
            sage: F1 = A3.free_graded_module([0])
            sage: g0 = F1.generator(0)
            sage: f = Hom(F0, F1)([A3.P(8)*g0, (A3.P(6,1))*g0])
            sage: f._resolve_kernel()
            Module morphism:
              From: Free graded left module on 5 generators over mod 3 Steenrod algebra, milnor basis
              To:   Free graded left module on 2 generators over mod 3 Steenrod algebra, milnor basis
              Defn: g[36] |--> P(1)*g[32]
                    g[44] |--> P(3)*g[32] + (2P(1))*g[40]
                    g[56] |--> P(6)*g[32] + P(0,1)*g[40]
                    g[64] |--> P(0,2)*g[32] + (2P(6))*g[40]
                    g[72] |--> P(6,1)*g[32]
        """
        return self._action(FPModuleMorphism._resolve_kernel, top_dim=top_dim, verbose=verbose)

    def _resolve_image(self, top_dim=None, verbose=False):
        r"""
        Resolve the image of this homomorphism by a free module.

        INPUT:

        - ``top_dim`` -- (optional) stop the computation at this degree; if
          not specified, this is determined using :meth:`profile`
        - ``verbose`` -- boolean (default: ``False``); whether log messages are
          printed

        OUTPUT: a homomorphism `j: F \rightarrow C` where `C` is the codomain
        of this homomorphism, `F` is free, and
        `\operatorname{im}(self) = \operatorname{im}(j)`.

        TESTS::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: F = SteenrodFPModule(A, [0,0])
            sage: L = SteenrodFPModule(A, [0,0], [[Sq(3),Sq(0,1)], [0,Sq(2)]])
            sage: f = Hom(F, L)([L([Sq(2),0]), L([0, Sq(2)])])
            sage: f._resolve_image()
            Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 2 generators and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[2] |--> Sq(2)*g[0, 0]
        """
        return self._action(FPModuleMorphism._resolve_image, top_dim=top_dim, verbose=verbose)

    def _action(self, method, *args, **kwds):
        r"""
        Changes the ground ring to a finite algebra, acts by the given method
        and changes back into the original ground ring before returning.

        TESTS::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: from sage.modules.fp_graded.morphism import FPModuleMorphism
            sage: A = SteenrodAlgebra(2)
            sage: F = SteenrodFPModule(A, [0])
            sage: L = SteenrodFPModule(A, [0], [[Sq(3)]])
            sage: f = Hom(F, L)([L([Sq(2)])])
            sage: f._action(FPModuleMorphism._resolve_image, verbose=True)
            Computing using the profile:
            (2, 1)
            Resolving the image in the range of dimensions [0, 8]: 0 1 2 3 4 5 6 7 8.
            Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[2] |--> Sq(2)*g[0]
        """
        small_profile = self.profile()

        if kwds.get('verbose', False):
            print('Computing using the profile:')
            print(small_profile)

        algebra = self.base_ring()

        # Choose a finite sub Hopf-algebra of the original algebra.
        finite_algebra = SteenrodAlgebra_generic(algebra.prime(), profile=small_profile)

        # Perform the chosen action on the module after having changed rings
        # to the finite algebra.
        fp_result = method(self.change_ring(finite_algebra), *args, **kwds)

        # Change back to the original algebra and also from FPModule
        # to SteenrodFPModule, and return the result.
        #
        # This is very clunky. Clean it up!
        f = fp_result.change_ring(self.base_ring())
        M = f.domain()
        N = f.codomain()
        new_values = [N.linear_combination(zip(N.generators(),
                                               v.dense_coefficient_list()))
                      for v in f.values()]
        return Hom(M, N)(new_values)


class SteenrodFreeModuleMorphism(FreeGradedModuleMorphism, SteenrodFPModuleMorphism):
    pass
