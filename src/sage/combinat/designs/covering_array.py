r"""
Covering Arrays (CA)

A Covering Array, denoted `CA(N;t,k,v)`, is an `N` by `k` array with
entries from a set of `v` elements with the property that in every
selection of `t` columns, every sequence of `t` elements appears in at
least one row.

An Orthogonal Array, denoted `OA(N;t,k,v)` is a covering array with the
property that every sequence of `t`-elements appears in exactly one row.
(See :mod:`sage.combinat.designs.orthogonal_arrays`).

This module collects methods relating to covering arrays, some of which
are inherited from orthogonal array methods. This module defines the
following functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.combinat.designs.designs_pyx.is_covering_array` | Check that an input list of lists is a `CA(N;t,k,v)`.
    :meth:`~sage.combinat.designs.covering_array.CA_relabel` | Return a relabelled version of the CA.
    :meth:`~sage.combinat.designs.covering_array.CA_standard_label` | Return a version of the CA relabelled to symbols `(0,\dots,n-1)`.

REFERENCES:

- [Colb2004]_

- [SMC2006]_

- [WC2007]_

AUTHORS:

- Aaron Dwyer and brett stevens (2022): initial version

"""

# **********************************************************************
#       Copyright (C) 2022 Aaron Dwyer <aarondwyer@cmail.carleton.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# **********************************************************************

from .orthogonal_arrays import OA_relabel, OA_standard_label
CA_relabel = OA_relabel
CA_standard_label = OA_standard_label
