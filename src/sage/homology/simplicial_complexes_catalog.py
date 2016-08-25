# -*- coding: utf-8 -*-
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

- :meth:`~sage.homology.examples.BarnetteSphere`
- :meth:`~sage.homology.examples.BrucknerGrunbaumSphere`
- :meth:`~sage.homology.examples.ChessboardComplex`
- :meth:`~sage.homology.examples.ComplexProjectivePlane`
- :meth:`~sage.homology.examples.DunceHat`
- :meth:`~sage.homology.examples.K3Surface`
- :meth:`~sage.homology.examples.KleinBottle`
- :meth:`~sage.homology.examples.MatchingComplex`
- :meth:`~sage.homology.examples.MooreSpace`
- :meth:`~sage.homology.examples.NotIConnectedGraphs`
- :meth:`~sage.homology.examples.PoincareHomologyThreeSphere`
- :meth:`~sage.homology.examples.PseudoQuaternionicProjectivePlane`
- :meth:`~sage.homology.examples.RandomComplex`
- :meth:`~sage.homology.examples.RandomTwoSphere`
- :meth:`~sage.homology.examples.RealProjectivePlane`
- :meth:`~sage.homology.examples.RealProjectiveSpace`
- :meth:`~sage.homology.examples.RudinBall`
- :meth:`~sage.homology.examples.ShiftedComplex`
- :meth:`~sage.homology.examples.Simplex`
- :meth:`~sage.homology.examples.Sphere`
- :meth:`~sage.homology.examples.SumComplex`
- :meth:`~sage.homology.examples.SurfaceOfGenus`
- :meth:`~sage.homology.examples.Torus`
- :meth:`~sage.homology.examples.ZieglerBall`

You can also get a list by typing ``simplicial_complexes.`` and hitting the
TAB key.

EXAMPLES::

    sage: S = simplicial_complexes.Sphere(2) # the 2-sphere
    sage: S.homology()
    {0: 0, 1: 0, 2: Z}
    sage: simplicial_complexes.SurfaceOfGenus(3)
    Triangulation of an orientable surface of genus 3
    sage: M4 = simplicial_complexes.MooreSpace(4)
    sage: M4.homology()
    {0: 0, 1: C4, 2: 0}
    sage: simplicial_complexes.MatchingComplex(6).homology()
    {0: 0, 1: Z^16, 2: 0}
"""

from sage.homology.examples import (Sphere, Simplex, Torus, ProjectivePlane,
        RealProjectivePlane, KleinBottle, SurfaceOfGenus, MooreSpace,
        ComplexProjectivePlane, PseudoQuaternionicProjectivePlane,
        PoincareHomologyThreeSphere, RealProjectiveSpace, K3Surface,
        BarnetteSphere, BrucknerGrunbaumSphere, NotIConnectedGraphs,
        MatchingComplex, ChessboardComplex, RandomComplex, SumComplex,
        RandomTwoSphere, ShiftedComplex, RudinBall, ZieglerBall, DunceHat)
