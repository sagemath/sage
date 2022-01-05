r"""
Catalog of simplicial sets

This provides pre-built simplicial sets:

- the `n`-sphere and `n`-dimensional real projective space, both (in
  theory) for any positive integer `n`. In practice, as `n` increases,
  it takes longer to construct these simplicial sets.

- the `n`-simplex and the horns obtained from it. As `n` increases, it
  takes *much* longer to construct these simplicial sets, because the
  number of nondegenerate simplices increases exponentially in `n`.
  For example, it is feasible to do
  ``simplicial_sets.RealProjectiveSpace(100)`` since it only has 101
  nondegenerate simplices, but ``simplicial_sets.Simplex(20)`` is
  probably a bad idea.

- `n`-dimensional complex projective space for `n \leq 4`

- the classifying space of a finite multiplicative group or monoid

- the torus and the Klein bottle

- the point

- the Hopf map: this is a pre-built morphism, from which one can
  extract its domain, codomain, mapping cone, etc.

All of these examples are accessible by typing
``simplicial_sets.NAME``, where ``NAME`` is the name of the
example. Type ``simplicial_sets.[TAB]`` for a complete list.

EXAMPLES::

    sage: RP10 = simplicial_sets.RealProjectiveSpace(8)
    sage: RP10.homology()
    {0: 0, 1: C2, 2: 0, 3: C2, 4: 0, 5: C2, 6: 0, 7: C2, 8: 0}

    sage: eta = simplicial_sets.HopfMap()
    sage: S3 = eta.domain()
    sage: S2 = eta.codomain()
    sage: S3.wedge(S2).homology()
    {0: 0, 1: 0, 2: Z, 3: Z}
"""

from .simplicial_set_examples import (Sphere, ClassifyingSpace,
                                      RealProjectiveSpace,
                                      KleinBottle, Torus,
                                      Simplex, Horn, Point,
                                      ComplexProjectiveSpace,
                                      HopfMap)
