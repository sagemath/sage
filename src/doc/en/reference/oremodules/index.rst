Modules over Ore rings
======================

Let `R` be a commutative ring, `\theta : K \to K` by a ring
endomorphism and `\partial : K \to K` be a `\theta`-derivation,
that is an additive map satisfying the following axiom

.. MATH::

    \partial(x y) = \theta(x) \partial(y) + \partial(x) y

The Ore polynomial ring associated to these data is 
`\mathcal S = R[X; \theta, \partial]`; its elements are the
usual polynomials over `R` but the multiplication is twisted
according to the rule

.. MATH::

    \partial(x y) = \theta(x) \partial(y) + \partial(x) y

We refer to :mod:`sage.rings.polynomial.ore_polynomial_ring.OrePolynomial`
for more details.

A Ore module over `(R, \theta, \partial)` is by definition a
module over `\mathcal S`; it is the same than a `R`-module `M`
equipped with an additive `f : M \to M` such that

.. MATH::

    f(a x) = \theta(a) f(x) + \partial(a) x

Such a map `f` is called a pseudomorphism
(see also :meth:`sage.modules.free_module.FreeModule_generic.pseudohom`).

SageMath provides support for creating and manipulating Ore
modules that are finite free over the base ring `R`.
This includes, in particular, Frobenius modules and modules
with connexions.

Modules, submodules and quotients
---------------------------------

.. toctree::
   :maxdepth: 1

   sage/modules/ore_module
   sage/modules/ore_module_element

Morphisms
---------

.. toctree::
   :maxdepth: 1

   sage/modules/ore_module_homspace
   sage/modules/ore_module_morphism
