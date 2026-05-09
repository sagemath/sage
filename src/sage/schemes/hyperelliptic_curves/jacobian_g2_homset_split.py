r"""
Rational point sets of Jacobians of genus-2 curves (split case)

AUTHORS:

- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#       Copyright (C) 2025 Sabrina Kunzweiler <sabrina.kunzweiler@math.u-bordeaux.fr>
#                     2025 Gareth Ma <grhkm21@gmail.com>
#                     2025 Giacomo Pope <giacomopope@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.hyperelliptic_curves.jacobian_homset_split import (
    HyperellipticJacobianHomsetSplit,
)


class HyperellipticJacobianHomsetSplit_g2(HyperellipticJacobianHomsetSplit):
    r"""
    Special class to handle optimisations for jacobian homset computations
    in genus two for hyperlliptic curves with a split model
    """

    pass
