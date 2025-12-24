r"""
Rational point sets of Jacobians of genus-2 curves (ramified case)

AUTHORS:

- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#       Copyright (C) 2025 Sabrina Kunzweiler, Gareth Ma, Giacomo Pope
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_ramified import (
    HyperellipticJacobianHomsetRamified,
)


class HyperellipticJacobianHomsetRamified_g2(HyperellipticJacobianHomsetRamified):
    """
    Special class to handle optimisations for jacobian homset computations
    in genus two for hyperlliptic curves with an ramified model
    """

    pass
