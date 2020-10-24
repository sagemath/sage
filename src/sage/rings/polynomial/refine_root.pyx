"""
Refine polynomial roots using Newton--Raphson

This is an implementation of the Newton--Raphson algorithm to
approximate roots of complex polynomials. The implementation
is based on interval arithmetic

AUTHORS:

- Carl Witty (2007-11-18): initial version
"""

#*****************************************************************************
#       Copyright (C) 2007 Carl Witty
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.real_mpfi import RealIntervalField
from sage.rings.complex_interval_field import ComplexIntervalField


def refine_root(ip, ipd, irt, fld):
    """
    We are given a polynomial and its derivative (with complex
    interval coefficients), an estimated root, and a complex interval
    field to use in computations.  We use interval arithmetic to
    refine the root and prove that we have in fact isolated a unique
    root.

    If we succeed, we return the isolated root; if we fail, we return
    None.

    EXAMPLES::

        sage: from sage.rings.polynomial.refine_root import refine_root
        sage: x = polygen(ZZ)
        sage: p = x^9 - 1
        sage: ip = CIF['x'](p); ip
        x^9 - 1
        sage: ipd = CIF['x'](p.derivative()); ipd
        9*x^8
        sage: irt = CIF(CC(cos(2*pi/9), sin(2*pi/9))); irt
        0.76604444311897802? + 0.64278760968653926?*I
        sage: ip(irt)
        0.?e-14 + 0.?e-14*I
        sage: ipd(irt)
        6.89439998807080? - 5.78508848717885?*I
        sage: refine_root(ip, ipd, irt, CIF)
        0.766044443118978? + 0.642787609686540?*I
    """

    # There has got to be a better way to do this, but I don't know
    # what it is...

    # We start with a basic fact: if we do an interval Newton-Raphson
    # step, and the refined interval is contained in the original interval,
    # then the refined interval contains exactly one root.

    # Unfortunately, our initial estimated root almost certainly does not
    # contain the actual root (our initial interval is a point, which
    # is exactly equal to whatever floating-point estimate we got from
    # the external solver).  So we need to do multiple Newton-Raphson
    # steps, and check this inclusion property on each step.

    # After a few steps of refinement, if the initial root estimate was
    # close to a root, we should have an essentially perfect interval
    # bound on the root (since Newton-Raphson has quadratic convergence),
    # unless either the real or imaginary component of the root is zero.
    # If the real or imaginary component is zero, then we could spend
    # a long time computing closer and closer approximations to that
    # component.  (This doesn't happen for non-zero components, because
    # of the imprecision of floating-point numbers combined with the
    # outward interval rounding; but close to zero, MPFI provides
    # extremely precise numbers.)

    # If the root is actually a real root, but we start with an imaginary
    # component, we can bounce back and forth between having a positive
    # and negative imaginary component, without ever hitting zero.
    # To deal with this, on every other Newton-Raphson step, instead of
    # replacing the old interval with the new one, we take the union.

    # If the containment check continues to fail many times in a row,
    # we give up and return None; we also return None if we detect
    # that the slope in our current interval is not bounded away
    # from zero at any step.

    # After every refinement step, we check to see if the real or
    # imaginary component of our interval includes zero.  If so, we
    # try setting it to exactly zero.  This gives us a good chance of
    # detecting real roots.  However, we do this replacement at most
    # once per component.

    refinement_steps = 10

    smashed_real = False
    smashed_imag = False

    for i in range(refinement_steps):
        slope = ipd(irt)
        if slope.contains_zero():
            return None
        center = fld(irt.center())
        val = ip(center)

        nirt = center - val / slope
        # print irt, nirt, (nirt in irt), nirt.diameter(), irt.diameter(), center, val, slope
        if nirt in irt and (nirt.diameter() >= irt.diameter() >> 3 or i >= 8):
            # If the new diameter is much less than the original diameter,
            # then we have not yet converged.  (Perhaps we were asked
            # for a particularly high-precision result.)  So we don't
            # return yet.
            return nirt

        if i & 1:
            irt = nirt
        else:
            irt = irt.union(nirt)
            # If we don't find a root after a while, try (approximately)
            # tripling the size of the region.
            if i >= 6:
                rD = irt.real().absolute_diameter()
                iD = irt.imag().absolute_diameter()
                md = max(rD, iD)
                md_intv = RealIntervalField(rD.prec())(-md, md)
                md_cintv = ComplexIntervalField(rD.prec())(md_intv, md_intv)
                irt = irt + md_cintv

        if not smashed_real and irt.real().contains_zero():
            irt = irt.parent()(0, irt.imag())
            smashed_real = True
        if not smashed_imag and irt.imag().contains_zero():
            irt = irt.parent()(irt.real(), 0)
            smashed_imag = True

    return None
