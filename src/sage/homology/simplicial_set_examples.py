# -*- coding: utf-8 -*-
r"""
Examples of simplicial sets: deprecated

The current version is :mod:`sage.topology.simplicial_set_examples`.
"""

from sage.misc.superseded import deprecated_function_alias
import sage.topology.simplicial_set_examples

for f in ['Nerve',
          'Sphere',
          'ClassifyingSpace',
          'RealProjectiveSpace',
          'KleinBottle',
          'Torus',
          'Simplex',
          'Empty',
          'Point',
          'Horn',
          'ComplexProjectiveSpace',
          'simplicial_data_from_kenzo_output',
          'HopfMap']:
    exec('{} = deprecated_function_alias(31925, sage.topology.simplicial_set_examples.{})'.format(f, f))
    
