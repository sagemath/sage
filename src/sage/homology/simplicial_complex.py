# -*- coding: utf-8 -*-
r"""
Finite simplicial complexes: deprecated

The current version is :mod:`sage.topology.simplicial_complexes`.
"""
from sage.misc.superseded import deprecated_function_alias
import sage.topology.simplicial_complex

for f in ['lattice_paths',
          'rename_vertex',
          'Simplex',
          'SimplicialComplex',
          'facets_for_RP4',
          'facets_for_K3']:
    exec('{} = deprecated_function_alias(31925, sage.topology.simplicial_complex.{})'.format(f, f))
