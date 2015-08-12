"""
Hyperelliptic curves of genus 2 over the rationals
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import hyperelliptic_g2_generic, hyperelliptic_rational_field

class HyperellipticCurve_g2_rational_field(
    hyperelliptic_g2_generic.HyperellipticCurve_g2_generic,
    hyperelliptic_rational_field.HyperellipticCurve_rational_field):
    pass
