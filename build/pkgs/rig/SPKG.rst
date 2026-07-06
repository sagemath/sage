rig: An optional GAP package
===================================

Description
-----------

Installing this SPKG will install the corresponding GAP package, but
before you can use them in Sage, they still have to be loaded into
either the GAP interface or libgap::

  sage: gap.eval('LoadPackage("rig")')  # optional - rig
  'true'
  sage: libgap.LoadPackage("srig")       # optional - rig
  true

Those correspond to::

  gap> LoadPackage("rig");

within the GAP interface and libgap, respectively.

Upstream Contact
----------------

See https://semigroups.github.io/Semigroups/

Dependencies
------------

-  GAP (a standard spkg) and  gap_packages (optional packages)

Notes
-----------
This is a GAP package for computations related to racks, quandles, knots, virtual knots, Nichols algebras.
(Author: L. Vendramin)
