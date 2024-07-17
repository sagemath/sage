.. -*- coding: utf-8 -*-
.. _steenrod_algebra_modules:

Steenrod Algebra Modules
========================

.. MODULEAUTHOR:: Robert R. Bruner, Michael J. Catanzaro, Sverre Lunoee--Nielsen, and Koen van Woerden

.. toctree::
   :maxdepth: 2

.. linkall

Let `p` be a prime number.  The mod `p` Steenrod algebra `A`
is a connected algebra over `\GF{p}`, the finite field of `p` elements.
All modules presented here will be defined over `A` or one of its sub-Hopf
algebras.  E.g.::

    sage: A = SteenrodAlgebra(p=2)

The constructor of the module class takes as arguments an ordered tuple of
degrees and the algebra over which the module is defined, together with an
optional set of relations::

    sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
    sage: F = SteenrodFPModule(A, [0, 1, 7]); F
    Free graded left module on 3 generators over mod 2 Steenrod algebra, milnor basis

Denote the module generators of an `A`-module `M` by `g_{d_1}, \ldots, g_{d_N}`,
where subscripts denote their degrees. A homogeneous relation of degree `n`
has the form

.. MATH::

    \sum_{i=1}^N a_i\cdot g_{d_i} = 0,

where the homogeneous coefficients `a_1, \ldots, a_N` lie in `A`, such that
`\deg(a_i) + \deg(g_{d_i}) = n` for `i = 1, \ldots, N`.  To create a module
with relations, the coefficients for each relation is given::

    sage: r1 = [Sq(8), Sq(7), 0]   # First relation
    sage: r2 = [Sq(7), 0, 1]       # Second relation
    sage: M = SteenrodFPModule(A, [0, 1, 7], relations=[r1, r2]); M
    Finitely presented left module on 3 generators and 2 relations
     over mod 2 Steenrod algebra, milnor basis

The resulting module will have three generators in the degrees we gave them::

    sage: M.generator_degrees()
    (0, 1, 7)

The connectivity of a module over a connected graded algebra is the minimum
degree of all its module generators.  Thus, if the module is non-trivial, the
connectivity is an integer::

    sage: M.connectivity()
    0

Each module is defined over a Steenrod algebra or some sub-Hopf algebra of it,
given by its base ring::

    sage: M.base_ring()
    mod 2 Steenrod algebra, milnor basis
    sage: SteenrodFPModule(SteenrodAlgebra(p=2,profile=(3,2,1)), [0]).base_ring()
    sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function
    [3, 2, 1]

.. NOTE::

    Calling ``algebra()`` will not return the desired algebra. Users should
    use the ``base_ring()`` method.

Module elements
---------------

Module elements are displayed in terms of generators, which by default
are called ``g[degree]``::

    sage: M.an_element(n=5)
    Sq(2,1)*g[0] + Sq(4)*g[1]

    sage: e = M.an_element(n=15); e
    Sq(0,0,0,1)*g[0] + Sq(1,2,1)*g[1] + Sq(8)*g[7]
    sage: e.dense_coefficient_list()
    [Sq(0,0,0,1), Sq(1,2,1), Sq(8)]

The generators are themselves elements of the module::

    sage: gens = M.generators(); gens
    (g[0], g[1], g[7])
    sage: gens[0] in M
    True

One can produce an element from a given set of algebra coefficients::

    sage: coeffs = [Sq(15), Sq(10)*Sq(1,1), Sq(8)]
    sage: x = M(coeffs); x
    Sq(15)*g[0] + (Sq(4,1,1)+Sq(7,0,1)+Sq(11,1))*g[1] + Sq(8)*g[7]

The module action produces new elements::

    sage: Sq(2) * x
    Sq(14,1)*g[0] + (Sq(7,1)+Sq(10))*g[7]

Each non-zero homogeneous element has a well-defined degree::

    sage: x.degree()
    15

but the zero element does not::

    sage: zero = M.zero(); zero
    0
    sage: zero.degree()
    Traceback (most recent call last):
    ...
    ValueError: the zero element does not have a well-defined degree

At this point it may be useful to point out that elements are not reduced to
a minimal representation when they are created.  A normalization can be
enforced, however::

    sage: g7 = M([0, 0, 1]); g7
    g[7]
    sage: g7.normalize()
    Sq(7)*g[0]
    sage: g7 == g7.normalize()
    True

    sage: m = M([Sq(7), 0, 0])
    sage: s = m + g7; s         # m and g7 are related by m = Sq(7)*g[0] = g[7],
    Sq(7)*g[0] + g[7]
    sage: s == 0                # so their sum should zero.
    True
    sage: s.normalize()         # Its normalized form is more revealing.
    0

For every integer `n`, the set of module elements of degree `n` form a
vector space over the ground field `\GF{p}`.  A basis for this vector space
can be computed::

    sage: M.basis_elements(7)
    (Sq(0,0,1)*g[0],
     Sq(1,2)*g[0],
     Sq(4,1)*g[0],
     Sq(7)*g[0],
     Sq(0,2)*g[1],
     Sq(3,1)*g[1],
     Sq(6)*g[1])

Note that the third generator `g_7` of degree 7 is apparently missing from the
basis above.  This is because of the relation `\operatorname{Sq}^7(g_0) = g_7`.

A vector space presentation can be produced::

    sage: M.vector_presentation(5)
    Vector space quotient V/W of dimension 4 over Finite Field of size 2 where
    V: Vector space of dimension 4 over Finite Field of size 2
    W: Vector space of degree 4 and dimension 0 over Finite Field of size 2
       Basis matrix:
       []

Given any element, its coordinates with respect to this basis can be computed::

    sage: x = M.an_element(7); x
    Sq(0,0,1)*g[0] + Sq(3,1)*g[1] + g[7]
    sage: v = x.vector_presentation(); v
    (1, 0, 0, 1, 0, 1, 0)

Going the other way, any element can be constructed by specifying its
coordinates::

    sage: x_ = M.element_from_coordinates((1, 0, 0, 1, 0, 1, 0), 7)
    sage: x_
    (Sq(0,0,1)+Sq(7))*g[0] + Sq(3,1)*g[1]
    sage: x_ == x
    True

Module homomorphisms
--------------------

Homomorphisms of `A`-modules `M\to N` are linear maps of their
underlying `\GF{p}`-vector spaces which commute with the `A`-module
structure. Homomorphisms are required to be homogeneneous but need not
be degree zero.

To create a homomorphism, first create the object modeling the set of all
such homomorphisms using the function ``Hom``::

    sage: Hko = SteenrodFPModule(A, [0], [[Sq(2)], [Sq(1)]])
    sage: homspace = Hom(Hko, Hko); homspace
    Set of Morphisms
     from Finitely presented left module on 1 generator and 2 relations
          over mod 2 Steenrod algebra, milnor basis
       to Finitely presented left module on 1 generator and 2 relations
          over mod 2 Steenrod algebra, milnor basis
       in Category of finitely presented graded modules
          over mod 2 Steenrod algebra, milnor basis

Just as with module elements, homomorphisms are created using the homspace.
The only argument is a list of elements in the codomain, giving the
images of the module generators of the domain::

    sage: gen = Hko.generator(0)  # the generator of the codomain module
    sage: values = [Sq(0, 0, 1) * gen]; values
    [Sq(0,0,1)*g[0]]
    sage: f = homspace(values)

The resulting homomorphism is the one sending the `i`-th generator of the
domain to the `i`-th codomain value given::

    sage: f
    Module endomorphism of Finitely presented left module on 1 generator
     and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> Sq(0,0,1)*g[0]

Homomorphisms can be evaluated on elements of the domain module::

    sage: v1 = f(Sq(4)*gen); v1
    Sq(4,0,1)*g[0]

    sage: v2 = f(Sq(2)*Sq(4)*gen); v2
    (Sq(3,1,1)+Sq(6,0,1))*g[0]

and they respect the module action::

    sage: f(Sq(4)*gen) == Sq(4)*f(gen)
    True

    sage: f(Sq(2)*Sq(4)*gen) == Sq(2)*Sq(4)*f(gen)
    True

Convenience methods exist for creating the trivial morphism::

    sage: x = Sq(4)*Sq(7)*gen
    sage: x == 0
    False
    sage: zero_map = homspace.zero(); zero_map
    Module endomorphism of Finitely presented left module on 1 generator
     and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> 0
    sage: zero_map(x)
    0
    sage: zero_map(x).is_zero()
    True

as well as the identity endomorphism::

    sage: one = Hom(Hko, Hko).identity(); one
    Module endomorphism of Finitely presented left module on 1 generator
     and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> g[0]
    sage: one.is_endomorphism()
    True
    sage: one(x) == x
    True
    sage: one.is_identity()
    True

Any non-trivial homomorphism has a well defined degree::

    sage: f.degree()
    7

but just as for module elements, the trivial homomorphism does not::

    sage: zero_map = homspace.zero()
    sage: zero_map.degree()
    Traceback (most recent call last):
    ...
    ValueError: the zero morphism does not have a well-defined degree

Any two homomorphisms can be added as long as they are of the same degree::

    sage: f1 = homspace([Hko([Sq(0,0,3) + Sq(0,2,0,1)])])
    sage: f2 = homspace([Hko([Sq(8,2,1)])])
    sage: (f1 + f2).is_zero()
    False
    sage: f1 + f2
    Module endomorphism of Finitely presented left module on 1 generator
     and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> (Sq(0,0,3)+Sq(0,2,0,1)+Sq(8,2,1))*g[0]

or when at least one of them is zero::

    sage: f + zero_map == f
    True

but not if they have different degrees::

    sage: F = SteenrodFPModule(A, [0])
    sage: b4 = Hom(F, F)([Sq(4) * F.generator(0)])
    sage: b8 = Hom(F, F)([Sq(8) * F.generator(0)])
    sage: b4 + b8
    Traceback (most recent call last):
    ...
    ValueError: morphisms do not have the same degree

Finally, additive inverses exist::

    sage: (f - f) == 0
    True

The restriction of a homomorphism to the vector space of `n`-dimensional module
elements is a linear transformation::

    sage: f_21 = f.vector_presentation(21); f_21
    Vector space morphism represented by the matrix:
    [1 0 0 0 0 0]
    [0 0 0 0 0 0]
    [1 0 0 0 0 0]
    Domain:   Vector space quotient V/W of dimension 3 over Finite Field of size 2 where
              V: Vector space of dimension 20 over Finite Field of size 2
              W: Vector space of degree 20 and dimension 17 over Finite Field of size 2
                 Basis matrix:
                 [1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
                 [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
                 [0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
                 [0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1]
                 [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1]
                 [0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 1 0 0 0 1]
                 [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1]
                 [0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 1]
                 [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
                 [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1]
                 [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
                 [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1]
                 [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1]
                 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1]
                 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1]
                 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
                 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1]
    Codomain: Vector space quotient V/W of dimension 6 over Finite Field of size 2 where
              V: Vector space of dimension 35 over Finite Field of size 2
              W: Vector space of degree 35 and dimension 29 over Finite Field of size 2
                 Basis matrix:
                 29 x 35 dense matrix over Finite Field of size 2

This is compatible with the vector presentations of its domain and codomain
modules::

    sage: f.domain() is Hko
    True
    sage: f.codomain() is Hko
    True
    sage: f_21.domain() is Hko.vector_presentation(21)
    True
    sage: f_21.codomain() is Hko.vector_presentation(21 + f.degree())
    True

Elements in the preimage of a homomorphism can be found::

    sage: f.solve(Sq(2)*Sq(4)*Sq(7)*gen)
    Sq(0,2)*g[0]

    sage: f.solve(Sq(8)*gen) is None
    True

Homomorphisms can be composed as expected::

    sage: g = homspace([Sq(0, 0, 0, 1)*gen]); g
    Module endomorphism of Finitely presented left module on 1 generator
     and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> Sq(0,0,0,1)*g[0]

    sage: g*f
    Module endomorphism of Finitely presented left module on 1 generator
     and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> Sq(0,0,1,1)*g[0]

    sage: one = homspace.identity()
    sage: f*one == f
    True


Homological algebra
^^^^^^^^^^^^^^^^^^^

The category of modules over `A` is Abelian, so kernels, images and
cokernels all exist and can be computed using :class:`the morphisms
<sage.modules.fp_graded.steenrod.morphism.SteenrodFPModuleMorphism>`.

.. NOTE::

    Methods for morphisms like

    - :meth:`~sage.modules.fp_graded.steenrod.morphism.SteenrodFPModuleMorphism.kernel_inclusion`,
    - :meth:`~sage.modules.fp_graded.steenrod.morphism.SteenrodFPModuleMorphism.cokernel_projection`,
    - :meth:`~sage.modules.fp_graded.steenrod.morphism.SteenrodFPModuleMorphism.image`,
    - :meth:`~sage.modules.fp_graded.morphism.FPModuleMorphism.homology`

    compute sub- and quotient modules related to homomorphisms, but
    they do not return instances of the module class.  Rather, they
    return the natural homomorphisms which connect these modules to
    the modules that gave rise to them.

    For example the function
    :meth:`~sage.modules.fp_graded.steenrod.morphism.SteenrodFPModuleMorphism.kernel_inclusion`
    returns an injective homomorphism which is onto the kernel
    submodule we asked it to compute, while the function
    :meth:`~sage.modules.fp_graded.steenrod.morphism.SteenrodFPModuleMorphism.cokernel_projection`
    provides a surjective homomorphism onto the cokernel module.

    In each case, getting a reference to the module instance requires calling
    :meth:`~sage.categories.morphism.Morphism.domain`
    or
    :meth:`~sage.categories.morphism.Morphism.codomain`
    on the returned homomorphism, depending on the case.

    Refer to each function's documentation for specific details.


Cokernels
^^^^^^^^^

In the following example, we define a cyclic module `H\ZZ` with one
relation in two ways: first explicitly, and then as the cokernel of a
homomorphism of free modules.  We then construct a candidate for an isomorphism
and check that it is both injective and surjective::

    sage: HZ = SteenrodFPModule(A, [0], [[Sq(1)]]); HZ
    Finitely presented left module on 1 generator and 1 relation
     over mod 2 Steenrod algebra, milnor basis

    sage: F = SteenrodFPModule(A, [0])
    sage: j = Hom(F, F)([Sq(1)*F.generator(0)])
    sage: coker = j.cokernel_projection()  # the natural quotient homomorphism onto the cokernel.
    sage: hz = coker.codomain(); hz
    Finitely presented left module on 1 generator and 1 relation
     over mod 2 Steenrod algebra, milnor basis

    sage: a = Hom(HZ, hz)([hz.generator(0)])
    sage: a.is_injective()
    True
    sage: a.is_surjective()
    True


Kernels
^^^^^^^

When computing the kernel of a homomorphism `f`, the result is an
injective homomorphism into the domain of `f`::

    sage: k = f.kernel_inclusion(); k
    Module morphism:
      From: Finitely presented left module on 1 generator and 3 relations
            over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 2 relations
            over mod 2 Steenrod algebra, milnor basis
      Defn: g[7] |--> Sq(0,0,1)*g[0]
    sage: k.codomain() == f.domain()
    True
    sage: k.is_injective()
    True

    sage: ker = k.domain()
    sage: ker
    Finitely presented left module on 1 generator and 3 relations
     over mod 2 Steenrod algebra, milnor basis

We can check that the injective image of `k` is the kernel of `f` by
showing that `f` factors as `h\circ c`, where `c` is the quotient map
to the cokernel of `k`, and `h` is injective::

    sage: K = k.codomain()        # We want to check that this really is the kernel of f.
    sage: coker = k.cokernel_projection()    # coker is the natural map: Hko -> coker(f) with kernel K.
    sage: h = Hom(coker.codomain(), Hko)(f.values())
    sage: h*coker == f            # Is K contained in ker(f) ?
    True
    sage: h.is_injective()        # Is ker(f) contained in K ?
    True


Images
^^^^^^

The method
:meth:`~sage.modules.fp_graded.steenrod.morphism.SteenrodFPModuleMorphism.image`
behaves similarly, returning an injective homomorphism with image
equal to the submodule `\operatorname{im}(f)`::

    sage: i = f.image(); i
    Module morphism:
      From: Finitely presented left module on 1 generator and 3 relations
            over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 2 relations
            over mod 2 Steenrod algebra, milnor basis
      Defn: g[7] |--> Sq(0,0,1)*g[0]
    sage: i.codomain() == f.codomain()
    True
    sage: i.is_injective()
    True

We can check that the injective image of `i` is the image of `f` by
lifting `f` over `i`, and showing that the lift is surjective::

    sage: f_ = f.lift(i); f_
    Module morphism:
      From: Finitely presented left module on 1 generator and 2 relations
            over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 3 relations
            over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> g[7]
    sage: i*f_ == f             # Is im(i) contained in im(f) ?
    True
    sage: f_.is_surjective()    # Is im(f) contained in im(i) ?
    True

When a pair of composable homomorphisms `g\circ f: M\to N\to L` satisfy the
condition `g\circ f = 0`, the sub-quotient `\ker(g) / \operatorname{im}(f)`
can be computed and is given by the natural quotient homomorphism with
domain `\ker(g)`::

    sage: f * f == 0        # Does the kernel of f contain the image of f ?
    True
    sage: K = f.kernel_inclusion()    # k: ker(f) -> Hko
    sage: h = f.homology(f) # h: ker(f) -> ker(f) / im(f)
    sage: h.codomain()      # This is the homology module.
    Finitely presented left module on 1 generator and 4 relations
     over mod 2 Steenrod algebra, milnor basis


Free resolutions
^^^^^^^^^^^^^^^^

Finally, free resolutions can be computed.  These calculations usually take
some time to complete, so it is usually a good idea to raise the verbose flag
to output progress information.

The following examples are taken from
`Michael Catanzaro's thesis <https://digitalcommons.wayne.edu/oa_theses/602/>`_
where the first version of this software appeared::

    sage: res = Hko.resolution(6, verbose=True)
    Computing f_1 (1/6)
    Computing f_2 (2/6)
    Computing using the profile:
    (2, 1)
    Resolving the kernel in the range of dimensions [1, 8]: 1 2 3 4 5 6 7 8.
    Computing f_3 (3/6)
    Computing using the profile:
    (2, 1)
    Resolving the kernel in the range of dimensions [2, 10]: 2 3 4 5 6 7 8 9 10.
    Computing f_4 (4/6)
    Computing using the profile:
    (2, 1)
    Resolving the kernel in the range of dimensions [3, 13]: 3 4 5 6 7 8 9 10 11 12 13.
    Computing f_5 (5/6)
    Computing using the profile:
    (2, 1)
    Resolving the kernel in the range of dimensions [4, 18]: 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18.
    Computing f_6 (6/6)
    Computing using the profile:
    (2, 1)
    Resolving the kernel in the range of dimensions [5, 20]: 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20.

The result of the calculation is a list of all the maps in the resolution::

    sage: [f.domain() for f in res]
    [Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 3 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 4 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 4 generators over mod 2 Steenrod algebra, milnor basis]

    sage: def is_complex(res):
    ....:     for i in range(len(res)-1):
    ....:         f = (res[i]*res[i+1])
    ....:         if not f.is_zero():
    ....:             return False
    ....:     return True
    ....:
    sage: is_complex(res)
    True

    sage: def is_exact(res):
    ....:     for i in range(len(res)-1):
    ....:         h = res[i].homology(res[i+1])
    ....:         if not h.codomain().is_trivial():
    ....:             return False
    ....:     return True

    sage: is_exact(res)
    True

    sage: [r.codomain().generator_degrees() for r in res]
    [(0,), (0,), (2, 1), (2, 4), (3, 7), (4, 8, 12), (5, 9, 13, 14)]

    sage: [r.values() for r in res]
    [(g[0],),
     (Sq(2)*g[0], Sq(1)*g[0]),
     (Sq(1)*g[1], Sq(3)*g[1] + Sq(2)*g[2]),
     (Sq(1)*g[2], Sq(2,1)*g[2] + Sq(3)*g[4]),
     (Sq(1)*g[3], Sq(2,1)*g[3] + Sq(1)*g[7], Sq(2,1)*g[7]),
     (Sq(1)*g[4],
      Sq(2,1)*g[4] + Sq(1)*g[8],
      Sq(2,1)*g[8] + Sq(1)*g[12],
      Sq(2)*g[12]),
     (Sq(1)*g[5],
      Sq(2,1)*g[5] + Sq(1)*g[9],
      Sq(2,1)*g[9] + Sq(1)*g[13],
      Sq(0,1)*g[13] + Sq(2)*g[14])]


Example: Counting lifts
-----------------------

In this more elaborate example we show how to find all possible lifts of a
particular homomorphism.  We will do this in two ways, and as a check of
validity, we will compare the results in the end.

We will work with the following modules over the mod 2
Steenrod algebra `A`:

.. MATH::
   :nowrap:

    \begin{align}
        H\ZZ &= A/A\cdot \operatorname{Sq}(1)\\
        Hko &= A/A\cdot \{\operatorname{Sq}(2), \operatorname{Sq}(1)\} \,.
    \end{align}

There is a natural projection `q: H\ZZ\to Hko`, and a non-trivial
endomorphism of degree 28, represented as a degree zero map
`f: \Sigma^{28}Hko\to Hko` that we define below.

The problem we will solve is to find all possible homomorphisms
`f': \Sigma^{28}Hko\to H\ZZ`, making the following diagram into a
commuting triangle:

.. MATH::

    H\ZZ\xrightarrow{q} Hko \xleftarrow{f} \Sigma^{28} Hko.

We begin by defining the modules and the homomorphisms `f` and `q`.
In the following, we let `L = \Sigma^{28}Hko`::

    sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
    sage: A = SteenrodAlgebra(2)
    sage: Hko = SteenrodFPModule(A, [0], [[Sq(2)],[Sq(1)]])
    sage: HZ = SteenrodFPModule(A, [0], [[Sq(1)]])
    sage: L = Hko.suspension(28)

The projection::

    sage: q = Hom(HZ, Hko)([Hko.generator(0)])
    sage: q
    Module morphism:
      From: Finitely presented left module on 1 generator and 1 relation
            over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 2 relations
            over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> g[0]

The map to lift over `q`::

    sage: f = Hom(L, Hko)([Sq(0,2,1,1)*Hko.generator(0)])
    sage: f
    Module morphism:
      From: Finitely presented left module on 1 generator and 2 relations
            over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 2 relations
            over mod 2 Steenrod algebra, milnor basis
      Defn: g[28] |--> Sq(0,2,1,1)*g[0]

    sage: f.is_zero()   # f is non-trivial.
    False

We will count the number of different lifts in two ways.  First, we will simply
compute the vector space of all possible maps `L \to H\ZZ`, and then check which
of those become `f` when composed with `q`::

    sage: basis = Hom(L, HZ).basis_elements(0)   # The basis for the vector space of degree 0 maps L -> HZ

    sage: from itertools import product
    sage: def from_coords(c):
    ....:     '''
    ....:     Create a linear combination of the three basis homomorphisms.
    ....:     '''
    ....:     return c[0]*basis[0] + c[1]*basis[1] + c[2]*basis[2]

    sage: for coords in product([0,1], repeat=3):
    ....:     print('%s: %s' % (coords, q*from_coords(coords) == f))
    (0, 0, 0): False
    (0, 0, 1): False
    (0, 1, 0): True
    (0, 1, 1): True
    (1, 0, 0): True
    (1, 0, 1): True
    (1, 1, 0): False
    (1, 1, 1): False

From this we conclude that four out of eight different homomorphisms
`L \to H\ZZ` are lifts of `f`::

    sage: lifts = [from_coords((0,1,0)),
    ....:          from_coords((0,1,1)),
    ....:          from_coords((1,0,0)),
    ....:          from_coords((1,0,1))]
    sage: lifts
    [Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> Sq(6,5,1)*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(6,5,1)+Sq(18,1,1))*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> Sq(10,1,0,1)*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(10,1,0,1)+Sq(18,1,1))*g[0]]

Alternatively we can use left-exactness of the functor
`\operatorname{Hom}_A(L, -)` to enumerate all possible lifts of `f`.
Start by finding a single lift of `f` over the projection `q`::

    sage: fl = f.lift(q); fl
    Module morphism:
      From: Finitely presented left module on 1 generator and 2 relations
            over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 1 relation
            over mod 2 Steenrod algebra, milnor basis
      Defn: g[28] |--> (Sq(4,3,0,1)+Sq(6,0,1,1)+Sq(7,2,0,1)+Sq(10,1,0,1))*g[0]

We verify that ``fl`` is indeed a lift::

    sage: q * fl == f
    True

There is an exact sequence

.. MATH::

    0 \to \operatorname{Hom}_A(L, \ker(q)) \xrightarrow{iK_*}
    \operatorname{Hom}_A(L, H\ZZ) \xrightarrow{q_*} \operatorname{Hom}_A(L, Hko),

which means that the indeterminacy of choosing a lift for
`f \in \operatorname{Hom}_A(L, Hko)` is represented by an element in
`\operatorname{Hom}_A(L,\ker(q))`.  Therefore, we can proceed to count the
number of lifts by computing this vector space of homomorphisms::

    sage: iK = q.kernel_inclusion()
    sage: K = iK.domain()
    sage: K.generator_degrees()
    (2,)
    sage: K.relations()
    (Sq(2)*g[2],)
    sage: ind = Hom(L, K).basis_elements(0); len(ind)
    2

So now we know that the vector space of indeterminacies is 2-dimensional over the
field of two elements.  This means that there are four distinct lifts of `f` over
`q`, and we can construct these by taking the one lift we already found, and add
to it all the different elements in the image of `iK_*`::

    sage: flift = [fl,
    ....:          fl + iK * ind[0],
    ....:          fl + iK * ind[1],
    ....:          fl + iK * (ind[0] + ind[1])]
    sage: flift
    [Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(4,3,0,1)+Sq(6,0,1,1)+Sq(7,2,0,1)+Sq(10,1,0,1))*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(0,7,1)+Sq(3,6,1)+Sq(4,1,3)+Sq(6,0,1,1)+Sq(6,5,1)+Sq(7,0,3))*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(4,3,0,1)+Sq(6,0,1,1)+Sq(7,2,0,1)+Sq(10,1,0,1)+Sq(12,3,1)+Sq(15,2,1)+Sq(18,1,1))*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(0,7,1)+Sq(3,6,1)+Sq(4,1,3)+Sq(6,0,1,1)+Sq(6,5,1)+Sq(7,0,3)+Sq(12,3,1)+Sq(15,2,1)+Sq(18,1,1))*g[0]]

As a test of correctness, we now compare the two sets of lifts.  As they stand,
it is not obvious that the lists ``flift`` and ``lifts`` are the same (up to a
re-ordering of list elements), so the following comparison is reassuring::

    sage: flift[0] == lifts[2]
    True
    sage: flift[1] == lifts[0]
    True
    sage: flift[2] == lifts[3]
    True
    sage: flift[3] == lifts[1]
    True
