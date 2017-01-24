# -*- coding: utf-8 -*-
r"""
Pointwise addition of dictionaries

This module is deprecated in favor of :mod:`sage.data_structures.blas_dict`.

EXAMPLES::

    sage: from sage.combinat.dict_addition import dict_addition, dict_linear_combination
    sage: D = { 0:1, 1:1 }; D
    {0: 1, 1: 1}

    sage: dict_addition( D for _ in range(5) )
    doctest:warning
    ...
    DeprecationWarning: dict_addition is deprecated. Please use sage.data_structures.blas_dict.sum instead.
    See http://trac.sagemath.org/20680 for details.
    {0: 5, 1: 5}

    sage: dict_linear_combination( (D,i) for i in range(5) )
    doctest:warning
    ...
    DeprecationWarning: dict_linear_combination is deprecated. Please use sage.data_structures.blas_dict.linear_combination instead.
    See http://trac.sagemath.org/20680 for details.
    {0: 10, 1: 10}
"""
#*****************************************************************************
#       Copyright (C) 2010 Christian Stump <christian.stump@univie.ac.at>
#                     2016 Travis Scrimshaw <tscrimsh@umn.edu>
#                     2016 Nicolas M. Thiéry <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecated_function_alias
import sage.data_structures.blas_dict as blas

dict_addition = deprecated_function_alias(20680, blas.sum)
dict_linear_combination = deprecated_function_alias(20680, blas.linear_combination)
