# sage.doctest: needs sage.graphs
"""
Catalog of simplicial complexes

There are two main types: manifolds and examples related to graph
theory.

For manifolds, there are functions defining the `n`-sphere for any
`n`, the torus, `n`-dimensional real projective space for any `n`, the
complex projective plane, surfaces of arbitrary genus, and some other
manifolds, all as simplicial complexes.

Aside from surfaces, this file also provides functions for
constructing some other simplicial complexes: the simplicial complex
of not-`i`-connected graphs on `n` vertices, the matching complex on n
vertices, the chessboard complex for an `n` by `i` chessboard, and
others.  These provide examples of large simplicial complexes; for
example, ``simplicial_complexes.NotIConnectedGraphs(7,2)`` has over a
million simplices.

All of these examples are accessible by typing
``simplicial_complexes.NAME``, where ``NAME`` is the name of the example.

- :meth:`~sage.topology.examples.BarnetteSphere`
- :meth:`~sage.topology.examples.BrucknerGrunbaumSphere`
- :meth:`~sage.topology.examples.ChessboardComplex`
- :meth:`~sage.topology.examples.ComplexProjectivePlane`
- :meth:`~sage.topology.examples.DunceHat`
- :meth:`~sage.topology.examples.FareyMap`
- :meth:`~sage.topology.examples.GenusSix`
- :meth:`~sage.topology.examples.K3Surface`
- :meth:`~sage.topology.examples.KleinBottle`
- :meth:`~sage.topology.examples.MatchingComplex`
- :meth:`~sage.topology.examples.MooreSpace`
- :meth:`~sage.topology.examples.NotIConnectedGraphs`
- :meth:`~sage.topology.examples.PoincareHomologyThreeSphere`
- :meth:`~sage.topology.examples.QuaternionicProjectivePlane`
- :meth:`~sage.topology.examples.RandomComplex`
- :meth:`~sage.topology.examples.RandomTwoSphere`
- :meth:`~sage.topology.examples.RealProjectivePlane`
- :meth:`~sage.topology.examples.RealProjectiveSpace`
- :meth:`~sage.topology.examples.RudinBall`
- :meth:`~sage.topology.examples.ShiftedComplex`
- :meth:`~sage.topology.examples.Simplex`
- :meth:`~sage.topology.examples.Sphere`
- :meth:`~sage.topology.examples.SumComplex`
- :meth:`~sage.topology.examples.SurfaceOfGenus`
- :meth:`~sage.topology.examples.Torus`
- :meth:`~sage.topology.examples.ZieglerBall`

You can also get a list by typing ``simplicial_complexes.`` and hitting the
:kbd:`Tab` key.

EXAMPLES::

    sage: S = simplicial_complexes.Sphere(2) # the 2-sphere
    sage: S.homology()                                                                  # needs sage.modules
    {0: 0, 1: 0, 2: Z}
    sage: simplicial_complexes.SurfaceOfGenus(3)
    Triangulation of an orientable surface of genus 3
    sage: M4 = simplicial_complexes.MooreSpace(4)
    sage: M4.homology()                                                                 # needs sage.modules
    {0: 0, 1: C4, 2: 0}
    sage: simplicial_complexes.MatchingComplex(6).homology()                            # needs sage.modules
    {0: 0, 1: Z^16, 2: 0}
"""

from sage.topology.simplicial_complex_examples import (Sphere, Simplex, Torus,
                                                       ProjectivePlane,
                                                       RealProjectivePlane, KleinBottle,
                                                       FareyMap, GenusSix,
                                                       SurfaceOfGenus,
                                                       MooreSpace,
                                                       ComplexProjectivePlane,
                                                       QuaternionicProjectivePlane,
                                                       PoincareHomologyThreeSphere,
                                                       RealProjectiveSpace, K3Surface,
                                                       BarnetteSphere,
                                                       BrucknerGrunbaumSphere,
                                                       NotIConnectedGraphs,
                                                       MatchingComplex,
                                                       ChessboardComplex, RandomComplex,
                                                       SumComplex,
                                                       RandomTwoSphere, ShiftedComplex,
                                                       RudinBall, ZieglerBall, DunceHat)
