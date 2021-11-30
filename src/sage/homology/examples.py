# -*- coding: utf-8 -*-
"""
Examples of simplicial complexes: deprecated

The current version is :mod:`sage.topology.simplicial_complex_examples`.
"""

from sage.misc.superseded import deprecated_function_alias
import sage.topology.simplicial_complex_examples

for f in ['facets_for_RP4',
          'facets_for_K3',
          'matching',
          'UniqueSimplicialComplex',
          'Sphere',
          'Simplex',
          'Torus',
          'RealProjectivePlane',
          'ProjectivePlane',
          'KleinBottle',
          'SurfaceOfGenus',
          'MooreSpace',
          'ComplexProjectivePlane',
          'PseudoQuaternionicProjectivePlane',
          'PoincareHomologyThreeSphere',
          'RealProjectiveSpace',
          'K3Surface',
          'BarnetteSphere',
          'BrucknerGrunbaumSphere',
          'NotIConnectedGraphs',
          'MatchingComplex',
          'ChessboardComplex',
          'RandomComplex',
          'SumComplex',
          'RandomTwoSphere',
          'ShiftedComplex',
          'RudinBall',
          'ZieglerBall',
          'DunceHat',
          'FareyMap']:
    exec('{} = deprecated_function_alias(31925, sage.topology.simplicial_complex_examples.{})'.format(f, f))
