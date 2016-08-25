"""
Hyperelliptic curves of genus 2 over a finite field
"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import hyperelliptic_finite_field, hyperelliptic_g2_generic

class HyperellipticCurve_g2_finite_field(
    hyperelliptic_g2_generic.HyperellipticCurve_g2_generic,
    hyperelliptic_finite_field.HyperellipticCurve_finite_field):
    pass

