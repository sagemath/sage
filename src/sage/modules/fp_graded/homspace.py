# sage.doctest: needs sage.combinat
r"""
Homsets of finitely presented graded modules

TESTS::

    sage: from sage.modules.fp_graded.module import FPModule
    sage: A = SteenrodAlgebra(2, profile=(3,2,1))
    sage: F = A.free_graded_module([1,3], names='c')
    sage: L = FPModule(A, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]], names='h')
    sage: homset = Hom(F, L); homset
    Set of Morphisms from Free graded left module on 2 generators ...
    sage: homset.an_element()
    Module morphism:
      From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      To:   Finitely presented left module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      Defn: c[1] |--> 0
            c[3] |--> Sq(1)*h[2]
    sage: f = homset([L((A.Sq(1), 1)), L((0, A.Sq(2)))]); f
    Module morphism:
      From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      To:   Finitely presented left module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      Defn: c[1] |--> Sq(1)*h[2] + h[3]
            c[3] |--> Sq(2)*h[3]
    sage: f.kernel_inclusion()
    Module morphism:
      From: Finitely presented left module on 2 generators and 1 relation over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      To:   Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      Defn: g[3] |--> c[3]
            g[4] |--> Sq(0,1)*c[1]

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

# *****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import Homset


class FPModuleHomspace(Homset):
    # FPModuleMorphism imports FPModuleHomspace, so this import should
    # not happen at the top level.
    from .morphism import FPModuleMorphism as Element

    def _element_constructor_(self, values):
        r"""
        Construct a morphism contained in ``self``.

        This function is not part of the public API, but is used by
        :meth:`Hom` method to create morphisms.

        INPUT:

        - ``values`` -- an iterable of FPElements of the codomain

        OUTPUT:

        A module homomorphism in this homspace sending the generators
        of the domain module to the given values.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FPModule(A2, [1,3])
            sage: L = FPModule(A2, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])

            sage: homset = Hom(F, L)
            sage: v1 = L([A2.Sq(1), 1])
            sage: v2 = L([0, A2.Sq(2)])
            sage: f = homset._element_constructor_([v1, v2])

        Rather than calling ``_element_constructor_`` explicitly, one
        should call it implicitly using the syntax::

            sage: f = homset([v1, v2]); f
            Module morphism:
              From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              To:   Finitely presented left module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[1] |--> Sq(1)*g[2] + g[3]
                    g[3] |--> Sq(2)*g[3]

        One can construct a homomorphism from another homomorphism::

            sage: g = homset(f)
            sage: f == g
            True

        And there is a convenient way of making the trivial homomorphism::

            sage: homset(0) == 0
            True

        TESTS::

            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = A2.free_graded_module((1,3))
            sage: L = A2.free_graded_module((2,5))
            sage: H = Hom(F, L)

            sage: values = (A2.Sq(4)*L.generator(0), A2.Sq(3)*L.generator(1))
            sage: f = H(values); f
            Module morphism:
              From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              To:   Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[1] |--> Sq(4)*g[2]
                    g[3] |--> Sq(3)*g[5]

            sage: H(0)
            Module morphism:
              From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              To:   Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[1] |--> 0
                    g[3] |--> 0
        """
        if isinstance(values, self.element_class):
            return values
        elif values == 0 or all(v.is_zero() for v in values):
            return self.zero()
        else:
            return self.element_class(self, values)

    def an_element(self, n=0):
        r"""
        Create a homomorphism belonging to ``self``.

        INPUT:

        - ``n`` -- (default: 0) an integer degree

        OUTPUT: a module homomorphism of degree ``n``

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: HZ = FPModule(A, [0], relations=[[Sq(1)]])

            sage: Hom(HZ, HZ).an_element(3)
            Module endomorphism of Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,1)*g[0]

        TESTS::

            sage: K = FPModule(A, [0, 0], [[Sq(2), 0]]) # Using a zero coefficient in the relations.
            sage: Hom(K, K).an_element(4)
            Module endomorphism of Finitely presented left module on 2 generators and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0, 0] |--> 0
                    g[0, 1] |--> Sq(4)*g[0, 0]

            sage: K = FPModule(A, [0, 0], [[Sq(2), 0], [0,0], [Sq(4), Sq(2)*Sq(2)]])
            sage: Hom(K, K).an_element(n=3)
            Module endomorphism of Finitely presented left module on 2 generators and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0, 0] |--> 0
                    g[0, 1] |--> Sq(0,1)*g[0, 0]

        An example involving free modules::

            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = A2.free_graded_module((1,3))
            sage: L = A2.free_graded_module((1,2))
            sage: H = Hom(F, L)
            sage: H.an_element()
            Module morphism:
              From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              To:   Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[1] |--> g[1]
                    g[3] |--> 0
        """
        return self._basis_elements(n, basis=False)

    def basis_elements(self, n):
        r"""
        Return a basis for the free module of degree ``n`` morphisms.

        INPUT:

        - ``n`` -- integer degree

        OUTPUT: a basis for the set of all module homomorphisms of degree ``n``

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: Hko = FPModule(A, [0], relations=[[Sq(2)], [Sq(1)]])

            sage: Hom(Hko, Hko).basis_elements(21)
            [Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
               Defn: g[0] |--> (Sq(0,0,3)+Sq(0,2,0,1))*g[0],
             Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
               Defn: g[0] |--> Sq(8,2,1)*g[0]]
        """
        return self._basis_elements(n, basis=True)

    def zero(self):
        r"""
        Create the trivial homomorphism in ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FPModule(A2, [1,3])
            sage: L = FPModule(A2, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])

            sage: z = Hom(F, L).zero()
            sage: z(F.an_element(5))
            0
            sage: z(F.an_element(23))
            0

        Example with free modules::

            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = A2.free_graded_module((1,3))
            sage: L = A2.free_graded_module((2,3))
            sage: H = Hom(F, L)
            sage: H.zero()
            Module morphism:
              From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              To:   Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[1] |--> 0
                    g[3] |--> 0
        """
        ngens = len(self.domain().generator_degrees())
        return self.element_class(self, [self.codomain().zero()] * ngens)

    def identity(self):
        r"""
        Return the identity homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: L = FPModule(A2, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])

            sage: one = Hom(L, L).identity(); one
            Module endomorphism of Finitely presented left module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[2] |--> g[2]
                    g[3] |--> g[3]

            sage: e = L.an_element(5)
            sage: e == one(e)
            True

        It is an error to call this function when the homset is not a
        set of endomorphisms::

            sage: F = FPModule(A2, [1,3])
            sage: Hom(F,L).identity()
            Traceback (most recent call last):
            ...
            TypeError: this homspace does not consist of endomorphisms

        An example with free graded modules::

            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: L = A2.free_graded_module((2,3))
            sage: H = Hom(L, L)
            sage: H.identity()
            Module endomorphism of Free graded left module on 2 generators
             over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
              Defn: g[2] |--> g[2]
                    g[3] |--> g[3]

        TESTS::

            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: L = A2.free_graded_module((2,3))
            sage: F = A2.free_graded_module((1,3))
            sage: H = Hom(F, L)
            sage: H.identity()
            Traceback (most recent call last):
            ...
            TypeError: this homspace does not consist of endomorphisms
        """
        if self.is_endomorphism_set():
            return self.element_class(self, self.codomain().generators())
        else:
            raise TypeError('this homspace does not consist of endomorphisms')

    def _basis_elements(self, n, basis):
        r"""
        Return a basis for the free module of degree ``n`` homomorphisms.

        This function is private for use by :meth:`basis_elements` and
        :meth:`an_element`.

        INPUT:

        - ``n`` -- integer degree
        - ``basis`` -- boolean; decide if a basis should be returned or just
          a single homomorphism

        OUTPUT:

        A basis for the set of all module homomorphisms of degree ``n``
        if ``basis`` is ``True``. Otherwise a single element is returned.
        In the latter case, this homomorphism is non-trivial if the vector
        space of all homomorphisms is non-trivial.

        TESTS::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: Hko = FPModule(A, [0], relations=[[Sq(2)], [Sq(1)]])
            sage: Hom(Hko, Hko)._basis_elements(21, basis=True)
            [Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
               Defn: g[0] |--> (Sq(0,0,3)+Sq(0,2,0,1))*g[0],
             Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
               Defn: g[0] |--> Sq(8,2,1)*g[0]]

            sage: Hom(Hko, Hko)._basis_elements(21, basis=False)
            Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> (Sq(0,0,3)+Sq(0,2,0,1))*g[0]

            sage: F = FPModule(A, [0])
            sage: Hom(F, Hko)._basis_elements(21, basis=False)
            Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,2,0,1)*g[0]

            sage: Hom(F, Hko)._basis_elements(21, basis=False)
            Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,2,0,1)*g[0]

        Test corner cases involving trivial modules::

            sage: F = FPModule(A, [0]) # A module without relations.
            sage: Z0 = FPModule(A, []) # A trivial module.
            sage: Z1 = FPModule(A, [0], [[1]]) # A trivial module with a redundant generator and relation.

            sage: Hom(FPModule(A, [-1]), F)._basis_elements(0, basis=True)
            []
            sage: Hom(FPModule(A, [-1]), F)._basis_elements(0, basis=False)
            Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[-1] |--> 0

            sage: from itertools import product
            sage: for D,C in product([(F, 'Free'), (Hko, 'Hko'), (Z0, 'Trivial'), (Z1, 'Trivial with redundant generator')], repeat=2):
            ....:     print('Hom(%s, %s):' % (D[1], C[1]))
            ....:     print('  basis==False:\n  %s' % Hom(D[0], C[0])._basis_elements(n=7, basis=False))
            ....:     print('  basis==True:\n  %s' % Hom(D[0], C[0])._basis_elements(n=7, basis=True))
            Hom(Free, Free):
              basis==False:
              Module endomorphism of Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,0,1)*g[0]
              basis==True:
              [Module endomorphism of Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,0,1)*g[0], Module endomorphism of Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(1,2)*g[0], Module endomorphism of Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(4,1)*g[0], Module endomorphism of Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(7)*g[0]]
            Hom(Free, Hko):
              basis==False:
              Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,0,1)*g[0]
              basis==True:
              [Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,0,1)*g[0]]
            Hom(Free, Trivial):
              basis==False:
              Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
            Hom(Free, Trivial with redundant generator):
              basis==False:
              Module morphism:
              From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
            Hom(Hko, Free):
              basis==False:
              Module morphism:
              From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
            Hom(Hko, Hko):
              basis==False:
              Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,0,1)*g[0]
              basis==True:
              [Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> Sq(0,0,1)*g[0]]
            Hom(Hko, Trivial):
              basis==False:
              Module morphism:
              From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
            Hom(Hko, Trivial with redundant generator):
              basis==False:
              Module morphism:
              From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
            Hom(Trivial, Free):
              basis==False:
              Module morphism:
              From: Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              basis==True:
              []
            Hom(Trivial, Hko):
              basis==False:
              Module morphism:
              From: Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              basis==True:
              []
            Hom(Trivial, Trivial):
              basis==False:
              Module endomorphism of Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
              basis==True:
              []
            Hom(Trivial, Trivial with redundant generator):
              basis==False:
              Module morphism:
              From: Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              basis==True:
              []
            Hom(Trivial with redundant generator, Free):
              basis==False:
              Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
            Hom(Trivial with redundant generator, Hko):
              basis==False:
              Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
            Hom(Trivial with redundant generator, Trivial):
              basis==False:
              Module morphism:
              From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              To:   Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
            Hom(Trivial with redundant generator, Trivial with redundant generator):
              basis==False:
              Module endomorphism of Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
              Defn: g[0] |--> 0
              basis==True:
              []
        """
        from .morphism import _create_relations_matrix

        M = self.domain()
        N = self.codomain()

        def _trivial_case():
            '''
            The return value if there are no non-trivial homomorphisms.
            '''
            if basis:
                # Since the free module of homomorphisms is trivial, the basis
                # is the empty set.
                return []
            else:
                # Since the free module of homomorphisms is trivial, it contains
                # only The trivial homomorphism
                return self.zero()

        # Deal with the trivial cases first.  Note that this covers the case
        # where the domain or codomain have no generators.
        if N.is_trivial() or M.is_trivial():
            return _trivial_case()

        # Then deal with the case where the domain has no relations.
        elif not M.has_relations():
            res = []
            num_generators = len(M.generators())
            for i, g in enumerate(M.generators()):
                # The i'th generator can go to any of these basis elements:
                base = N.basis_elements(g.degree() + n)
                for value in base:
                    values = [N.zero() if i != j else value for j in range(num_generators)]
                    res.append(self(values))
                    if not basis:
                        return res[0]

        else:
            # Note that this list is non-empty since we dealt with the trivial
            # case above.
            source_degs = [g.degree() + n for g in M.generators()]

            # Note that this list is non-empty since we dealt with the free
            # case above.
            target_degs = [r.degree() + n for r in M.relations()]

            block_matrix, R = _create_relations_matrix(
                N,
                [r.dense_coefficient_list() for r in M.relations()],
                source_degs, target_degs)

            ker = R.right_kernel()

            res = []
            for b in ker.basis():
                n = 0

                xs = []
                for j, X in enumerate(block_matrix[0]):
                    k = X.domain().dimension()
                    xs.append(N.element_from_coordinates(b[n:n + k], source_degs[j]))
                    n += k

                res.append(self(xs))
                if not basis:
                    return res[0]

        # If the code above found a non-trivial homomorphism and ``basis==False``,
        # it will have terminated by now.
        if not res:
            return _trivial_case()
        else:
            return res
