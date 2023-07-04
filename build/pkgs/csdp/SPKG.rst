csdp: Solver for semidefinite programs
======================================

Description
-----------

This is a fast SDP solver written in C, with a callable library namely,
an autotool'ed version of CSDP, by Brian Borchers, see
https://projects.coin-or.org/Csdp

License
-------

Common Public License Version 1.0


Upstream Contact
----------------

Dmitrii Pasechnik <dimpase+sage@gmail.com>

Special Update/Build Instructions
---------------------------------

csdp is an autotool'ed version of CSDP, see
https://projects.coin-or.org/Csdp, developed in its own repository at
https://github.com/dimpase/csdp.

To update to a new version, you need to bump the version number in
configure.ac and rerun autotools (autoreconf -fiv). Any changes should
be merged to the upstream repo.

The build is done with NOSHORTS variable defined; this makes it
compatible with packages, where NOSHORTS must be defined, e.g.
https://github.com/dimpase/pycsdp; also the Sage Cython interface needs
NOSHORTS defined.

Detailed steps to build the spkg are as follows. You need

-  git
-  autotools and libtool (the full autohell suite, version at least
   2.67)

With these ready:

-  ./spkg-src
-  copy the resulting csdp-<version>.tar.gz to SAGE_ROOT/upstream,
   or somewhere else appropriate
