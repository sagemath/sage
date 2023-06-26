r"""
Covering Arrays (CA)

A Covering Array, denoted CA`(N,k,v,t)`, is an `n` by `k` array with
entries from a set of `v` elements with the property that in every
selection of `t` columns, eachrow contains every sequence of
`t`-elements at least once.

An Orthogonal Array, denoted OA`(N,k,v,t)` is a covering array with the
property that each row contains every sequence of `t`-elements exactly
once

REFERENCES:

- [Colb2004]_

- [Sher2006]_

- [Wal2007]_

AUTHORS:

- Aaron Dwyer and brett stevens (2022): initial version

.. TODO::

    Implement various construction methods for CAs

Functions
---------

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
