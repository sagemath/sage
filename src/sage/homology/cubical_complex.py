# -*- coding: utf-8 -*-
r"""
Finite cubical complexes: deprecated

The current version is :mod:`sage.topology.cubical_complexes`.
"""

from sage.misc.superseded import deprecated_function_alias
import sage.topology.cubical_complex

Cube = deprecated_function_alias(31925, sage.topology.cubical_complex.Cube)
CubicalComplex = deprecated_function_alias(31925,
                  sage.topology.cubical_complex.CubicalComplex)
CubicalComplexExamples = deprecated_function_alias(31925,
                          sage.topology.cubical_complex.CubicalComplexExamples)
cubical_complexes = deprecated_function_alias(31925,
                     sage.topology.cubical_complex.cubical_complexes)
