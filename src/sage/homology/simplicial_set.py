# -*- coding: utf-8 -*-
r"""
Simplicial sets: deprecated

The current version is :mod:`sage.topology.simplicial_set`.
"""
from sage.misc.superseded import deprecated_function_alias
import sage.topology.simplicial_set

for f in ['AbstractSimplex_class',
          'NonDegenerateSimplex',
          'AbstractSimplex',
          'SimplicialSet_arbitrary',
          'SimplicialSet_finite',
          'SimplicialSet',
          'standardize_degeneracies',
          'all_degeneracies',
          'standardize_face_maps',
          'face_degeneracies',
          'shrink_simplicial_complex']:
    exec('{} = deprecated_function_alias(31925, sage.topology.simplicial_set.{})'.format(f, f))
