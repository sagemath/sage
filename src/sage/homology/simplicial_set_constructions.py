# -*- coding: utf-8 -*-
r"""
Methods of constructing simplicial sets: deprecated

The current version is :mod:`sage.topology.simplicial_set_constructions`.
"""
from sage.misc.superseded import deprecated_function_alias
import sage.topology.simplicial_set_constructions

for f in ['SubSimplicialSet',
          'PullbackOfSimplicialSets',
          'PullbackOfSimplicialSets_finite',
          'Factors',
          'ProductOfSimplicialSets',
          'ProductOfSimplicialSets_finite',
          'PushoutOfSimplicialSets',
          'PushoutOfSimplicialSets_finite',
          'QuotientOfSimplicialSet',
          'QuotientOfSimplicialSet_finite',
          'SmashProductOfSimplicialSets_finite',
          'WedgeOfSimplicialSets',
          'WedgeOfSimplicialSets_finite',
          'DisjointUnionOfSimplicialSets',
          'DisjointUnionOfSimplicialSets_finite',
          'ConeOfSimplicialSet',
          'ConeOfSimplicialSet_finite',
          'ReducedConeOfSimplicialSet',
          'ReducedConeOfSimplicialSet_finite',
          'SuspensionOfSimplicialSet',
          'SuspensionOfSimplicialSet_finite']:
    exec('{} = deprecated_function_alias(31925, sage.topology.simplicial_set_constructions.{})'.format(f, f))
