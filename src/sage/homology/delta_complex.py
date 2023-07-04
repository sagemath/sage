# -*- coding: utf-8 -*-
r"""
Finite Delta-complexes: deprecated

The current version is :mod:`sage.topology.delta_complexes`.
"""
from sage.misc.superseded import deprecated_function_alias
import sage.topology.delta_complex

DeltaComplex = deprecated_function_alias(31925,
                sage.topology.delta_complex.DeltaComplex)
DeltaComplexExamples = deprecated_function_alias(31925,
                        sage.topology.delta_complex.DeltaComplexExamples)
delta_complexes = deprecated_function_alias(31925,
                   sage.topology.delta_complex.delta_complexes)
