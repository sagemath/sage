r"""
Homomorphisms of finitely generated free graded left modules

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

# *****************************************************************************
#       Copyright (C) 2019 Robert R. Bruner <rrb@math.wayne.edu>
#                     and  Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.fp_graded.morphism import FPModuleMorphism


class FreeGradedModuleMorphism(FPModuleMorphism):
    r"""
    Create a homomorphism from a finitely generated free graded module
    to a graded module.

    INPUT:

    - ``parent`` -- a homspace in the category of finitely generated free
      modules

    - ``values`` -- list of elements in the codomain; each element
      corresponds (by their ordering) to a module generator in the domain

    EXAMPLES::

        sage: from sage.modules.fp_graded.free_module import FreeGradedModule
        sage: A = SteenrodAlgebra(2)
        sage: F1 = FreeGradedModule(A, (4,5), names='b')
        sage: F2 = FreeGradedModule(A, (3,4), names='c')
        sage: F3 = FreeGradedModule(A, (2,3), names='d')
        sage: H1 = Hom(F1, F2)
        sage: H2 = Hom(F2, F3)
        sage: f = H1((F2((Sq(4), 0)), F2((0, Sq(4)))))
        sage: g = H2((F3((Sq(2), 0)), F3((Sq(3), Sq(2)))))
        sage: g*f
        Module morphism:
          From: Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
          To:   Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
          Defn: b[4] |--> (Sq(0,2)+Sq(3,1)+Sq(6))*d[2]
                b[5] |--> (Sq(1,2)+Sq(7))*d[2] + (Sq(0,2)+Sq(3,1)+Sq(6))*d[3]

    TESTS:

    A non-example because the degree is not well-defined::

        sage: M = FreeGradedModule(A, (0, 0))
        sage: N = FreeGradedModule(A, (0,))
        sage: H = Hom(M, N)
        sage: g = N.generator(0)
        sage: H([Sq(1)*g, Sq(2)*g])
        Traceback (most recent call last):
        ...
        ValueError: ill-defined homomorphism: degrees do not match
    """

    def __init__(self, parent, values):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0, 0))
            sage: N = FreeGradedModule(A, (0,))
            sage: H = Hom(M, N)
            sage: g = N.generator(0)
            sage: TestSuite(H).run()
            sage: TestSuite(g).run()
        """
        from .free_homspace import FreeGradedModuleHomspace
        if not isinstance(parent, FreeGradedModuleHomspace):
            raise TypeError('the parent (%s) must be a f.p. free module homset' % parent)

        self._free_morphism = self
        FPModuleMorphism.__init__(self, parent, values, check=False)

        # Compute the degree.
        if all(v.is_zero() for v in self._values):
            # The zero homomorphism does not get a degree.
            degree = None
        else:
            degrees = []
            gen_deg = parent.domain().generator_degrees()
            for i, val in enumerate(self._values):
                if val:
                    x = val.degree()
                    xx = gen_deg[i]
                    degrees.append(x - xx)

            degree = min(degrees)
            if degree != max(degrees):
                raise ValueError('ill-defined homomorphism: degrees do not match')

        self._degree = degree

    def degree(self):
        r"""
        The degree of ``self``.

        OUTPUT:

        The degree of this homomorphism. Raise an error if this is
        the trivial homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: homspace = Hom(FreeGradedModule(A, (0,1)), FreeGradedModule(A, (0,)))
            sage: N = homspace.codomain()
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = homspace(values)
            sage: f.degree()
            5

        The zero homomorphism has no degree::

            sage: homspace.zero().degree()
            Traceback (most recent call last):
            ...
            ValueError: the zero morphism does not have a well-defined degree
        """
        if self._degree is None:
            # The zero morphism has no degree.
            raise ValueError("the zero morphism does not have a well-defined degree")
        return self._degree

    def __call__(self, x):
        r"""
        Evaluate the homomorphism at the given domain element ``x``.

        INPUT:

        - ``x`` -- an element of the domain of this morphism

        OUTPUT:

        The module element of the codomain which is the value of ``x``
        under this homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,1))
            sage: N = FreeGradedModule(A, (2,))
            sage: values = [Sq(5)*N.generator(0), Sq(3,1)*N.generator(0)]
            sage: f = Hom(M, N)(values)
            sage: f.__call__(M.generator(0))
            Sq(5)*g[2]
            sage: f.__call__(M.generator(1))
            Sq(3,1)*g[2]
        """
        if x.parent() != self.domain():
            raise ValueError('cannot evaluate morphism on element not in the domain')
        value = self.codomain().linear_combination(zip(self._values,
                                                       x.dense_coefficient_list()))
        return value

    def fp_module(self):
        r"""
        Create a finitely presented module from ``self``.

        OUTPUT: the finitely presented module with presentation equal to ``self``

        EXAMPLES::

            sage: A = SteenrodAlgebra(2)
            sage: F1 = A.free_graded_module([2])
            sage: F2 = A.free_graded_module([0])
            sage: v = F2([Sq(2)])
            sage: pres = Hom(F1, F2)([v])
            sage: M = pres.fp_module(); M
            Finitely presented left module on 1 generator and 1 relation over
             mod 2 Steenrod algebra, milnor basis
            sage: M.generator_degrees()
            (0,)
            sage: M.relations()
            (Sq(2)*g[0],)

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A = SteenrodAlgebra(2)
            sage: F1 = A.free_graded_module((2,))
            sage: F2 = FPModule(A, (0,), [[Sq(4)]])
            sage: v = F2([Sq(2)])
            sage: pres = Hom(F1, F2)([v])
            sage: pres.fp_module()
            Traceback (most recent call last):
            ...
            ValueError: this is not a morphism between free modules
        """
        if self.codomain().has_relations():
            raise ValueError("this is not a morphism between free modules")
        try:
            FPModule = self.base_ring()._fp_graded_module_class
        except AttributeError:
            from .module import FPModule
        return FPModule(self)
