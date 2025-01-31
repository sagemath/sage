# sage.doctest: needs sage.combinat
r"""
Finitely generated free graded left modules over connected graded algebras

Let `A` be a connected graded algebra. Some methods here require in
addition that `A` be an algebra over a field or a PID and that Sage
has a description of a basis for `A`.

For example, let `p` be a prime number. The mod `p` Steenrod algebra
`A_p` is a connected algebra over the finite field of `p` elements.
Many of the modules presented here will be defined over `A_p`, or one
of its sub-Hopf algebras.  E.g.::

    sage: A = SteenrodAlgebra(p=2)

However, the current implementation can use any connected graded
algebra that has a graded basis where each graded part is finite
dimensional. Another good family is the exterior algebras::

    sage: E.<x,y,z> = ExteriorAlgebra(QQ)

A free module is defined by the graded algebra and an ordered tuple
of degrees for the generators::

    sage: M = A.free_graded_module(generator_degrees=(0,1))
    sage: M
    Free graded left module on 2 generators over
     mod 2 Steenrod algebra, milnor basis

    sage: F.<a,b,c> = E.free_graded_module((0,3,6))
    sage: F
    Free graded left module on 3 generators over
     The exterior algebra of rank 3 over Rational Field

The resulting free modules will have generators in the degrees as specified::

    sage: M.generator_degrees()
    (0, 1)
    sage: F.generator_degrees()
    (0, 3, 6)

The default names for the generators are ``g[degree]`` if they are in
distinct degrees, ``g[degree, i]`` otherwise. They can be given other
names, as was done when creating the module ``F``::

    sage: M.generators()
    (g[0], g[1])
    sage: F.generators()
    (a, b, c)

The connectivity of a module over a connected graded algebra is the minimum
degree of all its module generators.  Thus, if the module is non-trivial, the
connectivity is an integer::

    sage: M.connectivity()
    0

Module elements
---------------

For an `A`-module with generators `\{g_i\}_{i=1}^N`, any homogeneous element
of degree `n` has the form

.. MATH::

    x = \sum_{i=1}^N a_i\cdot g_i\,,

where `a_i\in A_{n-\deg(g_i)}` for all `i`.  The ordered set `\{a_i\}`
is referred to as the coefficients of `x`.

You can produce module elements from a given set of coefficients::

    sage: coeffs = [Sq(5), Sq(1,1)]
    sage: x = M(coeffs); x
    Sq(5)*g[0] + Sq(1,1)*g[1]

You can also use the module action::

    sage: Sq(2) * x
    (Sq(4,1)+Sq(7))*g[0] + Sq(3,1)*g[1]

Each nonzero element has a well-defined degree::

    sage: x.degree()
    5

However the zero element does not::

    sage: zero = M.zero(); zero
    0
    sage: zero.degree()
    Traceback (most recent call last):
    ...
    ValueError: the zero element does not have a well-defined degree

Any two elements can be added as long as they are in the same degree::

    sage: y = M.an_element(5); y
    Sq(2,1)*g[0] + Sq(4)*g[1]
    sage: x + y
    (Sq(2,1)+Sq(5))*g[0] + (Sq(1,1)+Sq(4))*g[1]

or when at least one of them is zero::

    sage: x + zero == x
    True

Finally, additive inverses exist::

    sage: x - x
    0

For every integer `n`, the set of module elements of degree `n` form
a free module over the ground ring `k`.  A basis for this free module
can be computed::

    sage: M.basis_elements(5)
    (Sq(2,1)*g[0], Sq(5)*g[0], Sq(1,1)*g[1], Sq(4)*g[1])

together with a corresponding free module presentation::

    sage: M.vector_presentation(5)
    Vector space of dimension 4 over Finite Field of size 2

Given any element, its coordinates with respect to this basis can be computed::

    sage: v = x.vector_presentation(); v
    (0, 1, 1, 0)

Going the other way, any element can be constructed by specifying its
coordinates::

    sage: x_ = M.element_from_coordinates((0,1,1,0), 5)
    sage: x_
    Sq(5)*g[0] + Sq(1,1)*g[1]
    sage: x_ == x
    True


Module homomorphisms
--------------------

Homomorphisms of free graded `A`-modules `M\to N` are linear maps of their
underlying free `k`-module which commute with the `A`-module structure.

To create a homomorphism, first create the object modeling the set of all
such homomorphisms using the free function ``Hom``::

    sage: M = A.free_graded_module((0,1))
    sage: N.<c2> = A.free_graded_module((2,))
    sage: homspace = Hom(M, N); homspace
    Set of Morphisms from Free graded left module on 2 generators
      over mod 2 Steenrod algebra, milnor basis
     to Free graded left module on 1 generator
      over mod 2 Steenrod algebra, milnor basis
     in Category of finite dimensional graded modules with basis
      over mod 2 Steenrod algebra, milnor basis

Just as module elements, homomorphisms are created using the homspace
object. The only argument is a list of module elements in the codomain,
corresponding to the module generators of the domain::

    sage: values = [Sq(2)*c2, Sq(2)*Sq(1)*c2]
    sage: f = homspace(values)

The resulting homomorphism is the one sending the `i`-th generator of the
domain to the `i`-th codomain value given::

    sage: f
    Module morphism:
      From: Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
      To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> Sq(2)*c2
            g[1] |--> (Sq(0,1)+Sq(3))*c2

Convenience methods exist for creating the trivial morphism::

    sage: homspace.zero()
    Module morphism:
      From: Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
      To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> 0
            g[1] |--> 0

as well as the identity endomorphism::

    sage: Hom(M, M).identity()
    Module endomorphism of Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> g[0]
            g[1] |--> g[1]

Homomorphisms can be evaluated on elements of the domain module::

    sage: v1 = f(Sq(7)*M.generator(0)); v1
    Sq(3,2)*c2

    sage: v2 = f(Sq(17)*M.generator(1)); v2
    (Sq(11,3)+Sq(13,0,1)+Sq(17,1))*c2

and they respect the module action::

    sage: v1 == Sq(7)*f(M.generator(0))
    True

    sage: v2 == Sq(17)*f(M.generator(1))
    True

Any non-trivial homomorphism has a well-defined degree::

    sage: f.degree()
    4

but just as module elements, the trivial homomorphism does not::

    sage: zero_map = homspace.zero()
    sage: zero_map.degree()
    Traceback (most recent call last):
    ...
    ValueError: the zero morphism does not have a well-defined degree

Any two homomorphisms can be added as long as they are of the same degree::

    sage: f2 = homspace([Sq(2)*c2, Sq(3)*c2])
    sage: f + f2
    Module morphism:
      From: Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
      To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> 0
            g[1] |--> Sq(0,1)*c2

or when at least one of them is zero::

    sage: f + zero_map == f
    True

Finally, additive inverses exist::

    sage: f - f == 0
    True

The restriction of a homomorphism to the free module of `n`-dimensional
module elements is a linear transformation::

    sage: f_4 = f.vector_presentation(4); f_4
    Vector space morphism represented by the matrix:
    [0 1 0]
    [1 1 1]
    [0 1 0]
    [0 0 0]
    Domain: Vector space of dimension 4 over Finite Field of size 2
    Codomain: Vector space of dimension 3 over Finite Field of size 2

This is compatible with the vector presentations of its domain and codomain
modules::

    sage: f.domain() is M
    True
    sage: f.codomain() is N
    True
    sage: f_4.domain() is M.vector_presentation(4)
    True
    sage: f_4.codomain() is N.vector_presentation(4 + f.degree())
    True


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

from sage.misc.cachefunc import cached_method
from sage.modules.free_module import FreeModule
from sage.modules.fp_graded.free_element import FreeGradedModuleElement
from sage.rings.infinity import infinity
from sage.categories.graded_modules import GradedModules
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.homset import Hom
from sage.combinat.free_module import CombinatorialFreeModule


class FreeGradedModule(CombinatorialFreeModule):
    r"""
    Create a finitely generated free graded module over a connected
    graded algebra, with generators in specified degrees.

    INPUT:

    - ``algebra`` -- the graded connected algebra over which the module is
      defined; this algebra must be equipped with a graded basis

    - ``generator_degrees`` -- tuple of integers defining the number
      of generators of the module, and their degrees

    - ``names`` -- (optional) the names of the generators. If ``names``
      is a comma-separated string like ``'a, b, c'``, then those will
      be the names. Otherwise, for example if ``names`` is ``abc``,
      then the names will be ``abc(d,i)``.

    By default, if all generators are in distinct degrees, then the
    ``names`` of the generators will have the form ``g_{d}`` where
    ``d`` is the degree of the generator. If the degrees are not
    distinct, then the generators will be called ``g_{d,i}`` where
    ``d`` is the degree and ``i`` is its index in the list of
    generators in that degree.

    EXAMPLES::

        sage: from sage.modules.fp_graded.free_module import FreeGradedModule
        sage: E.<x,y,z> = ExteriorAlgebra(QQ)
        sage: M = FreeGradedModule(E, (-1,3))
        sage: M
        Free graded left module on 2 generators over
         The exterior algebra of rank 3 over Rational Field
        sage: M.generator_degrees()
        (-1, 3)
        sage: a, b = M.generators()
        sage: (x*y*b).degree()
        5

    ``names`` of generators::

        sage: M.generators()
        (g[-1], g[3])
        sage: FreeGradedModule(E, (0, 0, 2)).generators()
        (g[0, 0], g[0, 1], g[2, 0])
        sage: FreeGradedModule(E, (0, 0, 2), names='x, y, z').generators()
        (x, y, z)
        sage: FreeGradedModule(E, (0, 0, 2), names='xyz').generators()
        (xyz[0, 0], xyz[0, 1], xyz[2, 0])

    ``names`` can also be defined implicitly using Sage's ``M.<...>`` syntax::

        sage: A = SteenrodAlgebra(2)
        sage: M.<x,y,z> = FreeGradedModule(A, (-2,2,4))
        sage: M
        Free graded left module on 3 generators over
         mod 2 Steenrod algebra, milnor basis
        sage: M.gens()
        (x, y, z)
    """
    def __classcall__(cls, algebra, generator_degrees, category=None,
                      names=None, prefix=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: M1 = FreeGradedModule(E, [1, 0, 2], names='a,b,c')
            sage: M2.<a,b,c> = FreeGradedModule(E, (1, 0, 2))
            sage: M1 is M2
            True
        """
        if algebra.base_ring() not in PrincipalIdealDomains():
            raise ValueError('the ground ring of the algebra must be a PID')

        generator_degrees = tuple(generator_degrees)
        category = GradedModules(algebra).WithBasis().FiniteDimensional().or_subcategory(category)
        if names is not None:
            from sage.structure.category_object import normalize_names
            names = normalize_names(-1, names)
            if len(generator_degrees) > 1:
                if len(names) == 1:
                    if prefix is None:
                        prefix = names[0]
                    names = None  # if prefix is specified and takes priority
            if names is not None and len(names) != len(generator_degrees):
                raise ValueError("the names do not correspond to the generators")
        if prefix is None:
            prefix = 'g'
        return super().__classcall__(cls, algebra=algebra,
                                     generator_degrees=generator_degrees,
                                     category=category, names=names,
                                     prefix=prefix, **kwds)

    def __init__(self, algebra, generator_degrees, category, names=None, **kwds):
        r"""
        Create a finitely generated free graded module over a connected graded
        algebra.

        TESTS::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: TestSuite(FreeGradedModule(A, (-2,2,4))).run()
        """
        # If generator_degrees is [d_0, d_1, ...], then
        # the generators are indexed by (0,d_0), (1,d_1), ...
        keys = list(enumerate(generator_degrees))
        self._generator_degrees = generator_degrees

        keys = []
        degs_so_far = {}
        unique = True
        for i in generator_degrees:
            if i in degs_so_far:
                idx = degs_so_far[i] + 1
                degs_so_far[i] += 1
                unique = False
            else:
                idx = 0
                degs_so_far[i] = 0
            keys.append((i, idx))
        if unique:
            keys = [i[0] for i in keys]
        kwds['iterate_key'] = True

        # Call the base class constructor.
        CombinatorialFreeModule.__init__(self, algebra,
                                         basis_keys=keys,
                                         category=category,
                                         names=names,
                                         **kwds)

    Element = FreeGradedModuleElement

    def change_ring(self, algebra):
        r"""
        Change the base ring of ``self``.

        INPUT:

        - ``algebra`` -- a connected graded algebra

        OUTPUT:

        The free graded module over ``algebra`` defined with the same
        number of generators of the same degrees as ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))

            sage: M = FreeGradedModule(A, [0,1])
            sage: N = M.change_ring(A2); N
            Free graded left module on 2 generators over sub-Hopf algebra of
             mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]

        Changing back yields the original module::

            sage: N.change_ring(A) is M
            True
        """
        # We use the base class to avoid the category mixed one
        return type(self).__base__(algebra, self.generator_degrees(),
                                   prefix=self.prefix(), names=self._names)

    def _repr_(self):
        r"""
        Construct a string representation of ``self``.

        TESTS::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M
            Free graded left module on 3 generators over
             mod 2 Steenrod algebra, milnor basis
        """
        return ("Free graded left module on %s generator%s over %s"
                % (len(self._generator_degrees),
                   "" if len(self._generator_degrees) == 1 else "s",
                   self.base_ring()))

    def generator_degrees(self):
        r"""
        The degrees of the module generators.

        OUTPUT:

        A tuple containing the degrees of the generators for this
        module, in the order that the generators were given when
        ``self`` was constructed.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (-2,2,4))
            sage: M.generator_degrees()
            (-2, 2, 4)
        """
        return self._generator_degrees

    def is_trivial(self):
        r"""
        Return ``True`` if this module is trivial and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: FreeGradedModule(A, (-2,2,4)).is_trivial()
            False
            sage: FreeGradedModule(A, ()).is_trivial()
            True
        """
        return not self._generator_degrees

    def connectivity(self):
        r"""
        The connectivity of ``self``.

        OUTPUT:

        An integer equal to the minimal degree of all the generators, if
        this module is non-trivial.  Otherwise, `+\infty`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (-2,2,4))
            sage: M.connectivity()
            -2

        TESTS::

            sage: M = FreeGradedModule(A, ())
            sage: M.is_trivial()
            True
            sage: M.connectivity()
            +Infinity
        """
        return min(self.generator_degrees() + (infinity,))

    def _element_constructor_(self, coefficients):
        r"""
        Construct any element of ``self``.

        INPUT:

        - ``coefficients`` -- tuple of coefficient (i.e. elements of the
          algebra for this module), an element of FreeGradedModule, or the
          zero integer constant

        OUTPUT:

        An instance of the element class with coefficients from
        ``coefficients``, the element ``coefficients`` if it already
        was an element, or the zero module element.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M.<a,b,c> = FreeGradedModule(A, (0,2,4))

            sage: zero = M(0); zero
            0

            sage: e = M((Sq(4), Sq(2), 1)); e
            Sq(4)*a + Sq(2)*b + c

            sage: e is M(e)
            True
        """
        if isinstance(coefficients, self.element_class):
            return coefficients

        if not coefficients:
            return self.zero()

        A = self.base_ring()
        return self._from_dict({b: A(c) for (c, b) in zip(coefficients, self._indices) if c},
                               remove_zeros=False)

    def an_element(self, n=None):
        r"""
        Return an element of ``self``.

        This function chooses deterministically an element of the module
        in the given degree.

        INPUT:

        - ``n`` -- (optional) the degree of the element to construct

        OUTPUT:

        An element (of the given degree if specified).

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M.an_element(172)
            Sq(0,0,2,0,1,0,1)*g[0] + Sq(0,4,0,0,1,0,1)*g[2] + Sq(7,1,0,0,1,0,1)*g[4]

        Zero is the only element in the trivial module::

            sage: FreeGradedModule(A, ()).an_element()
            0
        """
        if not self._generator_degrees:
            return self.zero()

        if n is None:
            n = max(self.generator_degrees()) + 7

        coefficients = []
        for g in self.generator_degrees():
            basis = self.base_ring().basis(n - g) if n >= g else ()
            # All of the algebra generators in basis will bring the
            # module generator in dimension g to dimension
            # g + (topDimension - g) = topDimension.  Picking any one of them
            # will do, so we pick the one with index (g (mod l)).
            l = len(basis)
            if l:
                coefficients.append(basis[g % l])
            else:
                coefficients.append(self.base_ring().zero())

        return self(coefficients)

    # @cached_method
    def basis_elements(self, n):
        r"""
        Return a basis for the free module of degree ``n`` module elements.

        .. NOTE::

            Suppose ``self`` is a module over the graded algebra `A` with
            base ring `R`. This returns a basis as a free module over `R`,
            not a basis as a free module over `A`.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        A sequence of homogeneous module elements of degree ``n``, which
        is a basis for the free module of all degree ``n`` module elements.

        .. SEEALSO::

            :meth:`vector_presentation`, :meth:`element_from_coordinates`

        EXAMPLES::

            sage: A = SteenrodAlgebra(2)
            sage: M.<m0, m2, m4> = A.free_graded_module((0,2,4))
            sage: M.basis_elements(8)
            (Sq(1,0,1)*m0,
             Sq(2,2)*m0,
             Sq(5,1)*m0,
             Sq(8)*m0,
             Sq(0,2)*m2,
             Sq(3,1)*m2,
             Sq(6)*m2,
             Sq(1,1)*m4,
             Sq(4)*m4)
        """
        return tuple([self.term(self._indices[i], coeff)
                      for i in range(len(self._generator_degrees))
                      for coeff in self._basis_coeffs(n, i)])

    def _basis_coeffs(self, d, i):
        r"""
        Return a basis for the free module of degree ``d`` module elements
        corresponding to the ``i``-th generator.

        .. NOTE::

            Suppose ``self`` is a module over the graded algebra `A` with
            base ring `R`. This returns a basis as a free module over `R`,
            not a basis as a free module over `A`.

        INPUT:

        - ``d`` -- integer; the degree
        - ``i`` -- integer; the factor

        EXAMPLES::

            sage: A = SteenrodAlgebra(2)
            sage: M = A.free_graded_module((0,1))
            sage: M._basis_coeffs(3, 0)
            (Sq(0,1), Sq(3))
            sage: M._basis_coeffs(3, 1)
            (Sq(2),)
        """
        return self._cached_basis_coeffs(d - self._generator_degrees[i])

    @cached_method
    def _cached_basis_coeffs(self, d):
        """
        Return the basis for the degree ``d`` part of the base algebra.

        EXAMPLES::

            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: M = E.free_graded_module((0,0,1))
            sage: M._cached_basis_coeffs(0)
            (1,)
            sage: M._cached_basis_coeffs(1)
            (x, y)
            sage: M._cached_basis_coeffs(2)
            (x*y,)
            sage: M._cached_basis_coeffs(3)
            ()
        """
        return tuple(self.base_ring().basis(d))

    @cached_method
    def element_from_coordinates(self, coordinates, n):
        r"""
        The module element of degree ``n`` having the given coordinates
        with respect to the basis of module elements given by
        :meth:`basis_elements`.

        INPUT:

        - ``coordinates`` -- a sequence of elements of the ground ring
        - ``n`` -- integer

        OUTPUT: a module element of degree ``n``

        .. SEEALSO::

            :meth:`vector_presentation`, and :meth:`basis_elements`.

        EXAMPLES::

            sage: A = SteenrodAlgebra(2)
            sage: M = A.free_graded_module((0,1))
            sage: x = M.element_from_coordinates((0,1,0,1), 5); x
            Sq(5)*g[0] + Sq(4)*g[1]
            sage: basis = M.basis_elements(5)
            sage: y = 0*basis[0] + 1*basis[1] + 0*basis[2] + 1*basis[3]
            sage: x == y
            True

            sage: M.element_from_coordinates((0,0,0,0), 5)
            0
        """
        D = self.vector_presentation(n).dimension()
        if len(coordinates) != D:
            raise ValueError('the given coordinate vector has incorrect length (%d); '
                             'it should have length %d' % (len(coordinates), D))

        # Performance testing using this real life example:
        #
        # sage: A = SteenrodAlgebra(2)
        # sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
        # sage: rels = [[Sq(1),0,0,0], [Sq(2),0,0,0], [Sq(4),0,0,0], [Sq(8),0,0,0],
        # ....:         [0,Sq(1),0,0], [0,Sq(2),0,0], [0,Sq(4),0,0], [Sq(31),Sq(14),0,0],
        # ....:         [0,Sq(20),0,0], [0,0,Sq(1),0], [0,0,Sq(2),0], [0,Sq(31),Sq(6),0],
        # ....:         [0,0,Sq(8),0], [0,0,0,Sq(1)], [0,0,Sq(31),Sq(2)], [0,0,0,Sq(4)], [0,0,0,Sq(8)] ]
        # sage: M = SteenrodFPModule(A, [0, 17, 42, 71], relations=rels)
        # sage: res = M.resolution(2, top_dim=30, verbose=True)
        #
        # This function was called a total of 2897 times during the computation,
        # and the total running time of the entire computation dropped from
        # 57 to 21 seconds by adding the optimization.

        ret = {}
        A = self.base_ring()
        j = 0
        for i, key in enumerate(self._indices):
            B = self._basis_coeffs(n, i)
            coeff = A.linear_combination((b, coordinates[j + ind])
                                         for ind, b in enumerate(B))
            if coeff:
                ret[key] = coeff
            j += len(B)

        if not ret:
            return self.zero()
        return self.element_class(self, ret)

    @cached_method
    def vector_presentation(self, n):
        r"""
        Return a free module over the ground ring of the module algebra
        isomorphic to the degree ``n`` elements of ``self``.

        Let `\mathcal{k}` be the ground ring of the algebra over this module
        is defined, and let `M_n` be the free module of module elements of
        degree ``n``.

        The return value of this function is the free module
        `\mathcal{k}^{r}` where `r = dim(M_n)`.

        The isomorphism between `k^{r}` and `M_n` is given by the
        bijection taking the standard basis element `e_i` to the `i`-th
        element of the array returned by :meth:`basis_elements`.

        INPUT:

        - ``n`` -- integer degree

        OUTPUT:

        A free module over the ground ring of the algebra over which
        ``self`` is defined, isomorphic to the free module of module
        elements of degree ``n``.

        .. SEEALSO::

            :meth:`basis_elements`, :meth:`element_from_coordinates`

        EXAMPLES::

            sage: A1 = SteenrodAlgebra(2, profile=[2,1])
            sage: M.<x> = A1.free_graded_module((0,))
            sage: M.vector_presentation(3)
            Vector space of dimension 2 over Finite Field of size 2
            sage: M.basis_elements(3)
            (Sq(0,1)*x, Sq(3)*x)
            sage: [M.vector_presentation(i).dimension() for i in range(-2, 9)]
            [0, 0, 1, 1, 1, 2, 1, 1, 1, 0, 0]

        TESTS::

            sage: A = SteenrodAlgebra(2)
            sage: M = A.free_graded_module((0,2,4))
            sage: V = M[4]; V
            Vector space of dimension 4 over Finite Field of size 2
            sage: V.dimension()
            4
        """
        m = len(self._generator_degrees)
        return FreeModule(self.base_ring().base_ring(), sum(len(self._basis_coeffs(n, i))
                                                            for i in range(m)))

    __getitem__ = vector_presentation

    def generator(self, index):
        r"""
        Return the module generator with the given index.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M.generator(0)
            g[0]
            sage: M.generator(1)
            g[2]
            sage: M.generator(2)
            g[4]
        """
        try:
            return self.gens()[index]
        except IndexError:
            raise ValueError('the parent module has generators in the index '
                             'range [0, %s]; generator %s does not exist' %
                             (len(self.generator_degrees()) - 1, index))

    gen = generator

    def generators(self):
        r"""
        Return all the module generators.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (-2,1))
            sage: M.generators()
            (g[-2], g[1])
        """
        return self.gens()

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
        from .free_homspace import FreeGradedModuleHomspace
        return FreeGradedModuleHomspace(self, Y, category)

    def suspension(self, t):
        r"""
        Suspend ``self`` by the given degree ``t``.

        INPUT:

        - ``t`` -- integer

        OUTPUT:

        A module which is isomorphic to this module by a shift
        of degrees by the integer ``t``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeGradedModule(A, (0,2,4))
            sage: M.suspension(4).generator_degrees()
            (4, 6, 8)
            sage: M.suspension(-4).generator_degrees()
            (-4, -2, 0)
        """
        degs = tuple(g + t for g in self.generator_degrees())
        return FreeGradedModule(algebra=self.base_ring(),
                                generator_degrees=degs)

    def has_relations(self):
        r"""
        Return ``False`` as this has no relations.

        This is for compatibility with
        :class:`~sage.modules.fp_graded.module.FPModule`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: F = FreeGradedModule(A, (-2,2,4))
            sage: F.has_relations()
            False
        """
        return False

    def relations(self):
        r"""
        Return the relations of ``self``, which is ``()``.

        This is for compatibility with
        :class:`~sage.modules.fp_graded.module.FPModule`.

        EXAMPLES::

            sage: from sage.modules.fp_graded.free_module import FreeGradedModule
            sage: A = SteenrodAlgebra(2)
            sage: F = FreeGradedModule(A, (-2,2,4))
            sage: F.relations()
            ()
        """
        return ()

    def resolution(self, k, top_dim=None, verbose=False):
        r"""
        Return a free resolution of ``self`` of length ``k``.

        Since ``self`` is free, the initial map in the resolution will
        be the identity, and the rest of the maps will be zero.

        INPUT:

        - ``k`` -- nonnegative integer
        - ``top_dim`` -- stop the computation at this degree. Ignored,
          for compatibility with
          :meth:`sage.modules.fp_graded.module.FPModule.resolution`.
        - ``verbose`` -- boolean (default: ``False``); controls whether log
          messages should be emitted

        OUTPUT:

        A list of homomorphisms `[1_M, 0, 0, \ldots, 0]` consisting of
        the identity map on this module followed by zero maps. Other
        than this module, the other modules in the resolution will be
        zero.

        EXAMPLES::

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: M = E.free_graded_module((1,2))
            sage: M.resolution(0)
            [Module endomorphism of Free graded left module on 2 generators over The exterior algebra of rank 3 over Rational Field
               Defn: g[1] |--> g[1]
                     g[2] |--> g[2]]
            sage: M.resolution(1)
            [Module endomorphism of Free graded left module on 2 generators over The exterior algebra of rank 3 over Rational Field
               Defn: g[1] |--> g[1]
                     g[2] |--> g[2],
             Module morphism:
               From: Free graded left module on 0 generators over The exterior algebra of rank 3 over Rational Field
               To:   Free graded left module on 2 generators over The exterior algebra of rank 3 over Rational Field]
            sage: M.resolution(4)
            [Module endomorphism of Free graded left module on 2 generators over The exterior algebra of rank 3 over Rational Field
               Defn: g[1] |--> g[1]
                     g[2] |--> g[2],
             Module morphism:
               From: Free graded left module on 0 generators over The exterior algebra of rank 3 over Rational Field
               To:   Free graded left module on 2 generators over The exterior algebra of rank 3 over Rational Field,
             Module endomorphism of Free graded left module on 0 generators over The exterior algebra of rank 3 over Rational Field,
             Module endomorphism of Free graded left module on 0 generators over The exterior algebra of rank 3 over Rational Field,
             Module endomorphism of Free graded left module on 0 generators over The exterior algebra of rank 3 over Rational Field]
        """
        if k < 0:
            raise ValueError('the length of the resolution must be nonnegative')

        # The first map \epsilon is the identity map
        ret_complex = [Hom(self, self).identity()]

        if k == 0:
            return ret_complex

        # All subsequent maps are trivial since self is free
        T = self.base_ring().free_graded_module(())
        ret_complex.append(Hom(T, self).zero())

        if k == 1:
            return ret_complex

        return ret_complex + [Hom(T, T).zero()] * (k - 1)

    def minimal_presentation(self, top_dim=None, verbose=False):
        r"""
        Return a minimal presentation of ``self``.

        OUTPUT: the identity morphism as ``self`` is free

        EXAMPLES::

            sage: from sage.modules.fp_graded.module import FPModule
            sage: A2 = SteenrodAlgebra(2)

            sage: M = A2.free_graded_module([0,1])
            sage: M.minimal_presentation().is_identity()
            True
        """
        return Hom(self, self).identity()
