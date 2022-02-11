zn_poly: C library for polynomial arithmetic in Z/nZ[x]
=======================================================

Description
-----------

zn_poly is a C library for polynomial arithmetic in Z/nZ[x], where n is
any modulus that fits into an unsigned long.

Website: https://gitlab.com/sagemath/zn_poly

Note: Original website is at https://web.maths.unsw.edu.au/~davidharvey/code/zn_poly/ but is
no longer maintained. Sage maintains an "official" continuation of the
project at the above link.

License
-------

GPL V2 or V3. Some of the code has been copied from other projects - see
the file src/COPYING for details.


Upstream Contact
----------------

-  David Harvey
-  \E. M. Bray <erik.m.bray@gmail.com>

Dependencies
------------

-  GMP/MPIR
-  (some) Python (to create the Makefile)
-  GNU patch
-  NTL apparently only if we configured zn_poly differently (same for
   FLINT)


Special Update/Build Instructions
---------------------------------

-  Make sure the patches still apply.

   Especially changes in ``makemakefile.py`` may also require changes to
   ``spkg-install`` (and perhaps also ``spkg-check``).

-  There's also a ``--use-flint`` option to ``configure``; no idea what
   it does,
   and we currently don't use it either.

-  TODO:
-  Use ``make install`` instead of manually "installing" (copying and
   symlinking) the [shared] libraries and header files. This requires
   further
   tweaking of ``makemakefile.py``, since it currently only installs a
   static
   library and the headers.

-  If everything's fine, i.e., no problems arise, some comments and
   especially some code I currently just commented out can certainly be removed.
   (-leif, 04/2012)

-  The version number "0.9.p11" is used as a doctest in the function
   package_versions in sage/misc/packages.py, so if this package gets
   upgraded, that doctest needs to be changed.

Patches
~~~~~~~

-  All patches from Sage have been merged into upstream. These include:
-  makemakefile.py.patch:

   Improves the Python script creating the Makeefile for better use at
   least within Sage; see patch for details. (Last modified at #12433,
   which added and changed a lot.)

-  profiler.c.patch, zn_poly.h.patch:

   Fix potential redefinition of ``ulong`` (in combination with other
   headers).

-  mpn_mulmid-tune.c.patch, mulmid-tune.c.patch, mul-tune.c.patch:

   Fix "jump into scope of identifier with variably modified type"
   errors. (See #8771).

-  mpn_mulmid-test.c.patch:

   Fix a potential problem when the value of ZNP_mpn_smp_kara_thresh is
   SIZE_MAX, this is usually irrealistic but can happen at least on
   linux on power7 with gcc-4.7.1 (see #14098).

-  fix_fudge_factor_in_nuss-test.c.patch:

   As the name says; fix provided by upstream (David Harvey); see
   #13947.
