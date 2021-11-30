gf2x: Fast arithmetic in GF(2)[x] and searching for irreducible/primitive trinomials
====================================================================================

Description
-----------

gf2x is a C/C++ software package containing routines for fast arithmetic
in GF(2)[x] (multiplication, squaring, GCD) and searching for
irreducible/primitive trinomials.

Website: http://gf2x.gforge.inria.fr/

License
-------

-  GNU GPLv2+.


Upstream Contact
----------------

-  Richard Brent
-  Pierrick Gaudry
-  Emmanuel Thomé
-  Paul Zimmermann

Dependencies
------------

-  None


Special Update/Build Instructions
---------------------------------

-  As some patches touch config/acinclude.m4, we have to touch
   aclocal.m4,
   configure, Makefile.in and gf2x/gf2x-config.h.in to prevent autotools
   to try to regenerate these files.

Patches
~~~~~~~

-  0001-Trac-15014-Let-gf2x-build-a-shared-library-on-Cygwin.patch: pass
   -no-undefined flag to libtool.
-  0002-tr-portability.patch: backport upstream fix for non-portable tr
   use
-  0003-Improve-detection-of-sse2-support.patch: backport upstream
   improved check for sse2

-  0004-Add-disable-hardware-specific-code.patch: add option
   -disable-hardware-specific-code to build system. This is partly
   backported from upstream.

-  0005-Update-autotooled-files.patch: the above patches make changes to
   code used by autotools for generation of the build system. This
   patches
   those files, so that autotools need not be installed.

-  0006-Fix_make_check_not_failing_on_errors.patch: (upstream patch)
   Fix bug in shell script such that 'make check' always fails upon
   errors.
