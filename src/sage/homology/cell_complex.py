# -*- coding: utf-8 -*-
r"""
Generic cell complexes: deprecated

The current version is :mod:`sage.topology.cell_complexes`.
"""

from sage.misc.superseded import deprecated_function_alias
import sage.topology.cell_complex

GenericCellComplex = deprecated_function_alias(31925,
                      sage.topology.cell_complex.GenericCellComplex)
