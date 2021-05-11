mpir: Multiple precision integers and rationals (fork of GMP)
=============================================================

Description
-----------

Multiple Precision Integers and Rationals.

MPIR is an open source multiprecision integer library derived from
version 5.0.1 of the GMP (GNU Multi Precision) project (which was
licensed LGPL v2+).

See http://www.mpir.org

License
-------

-  LGPL V3+


Upstream Contact
----------------

-  The Google group mpir-devel
-  thempirteam@googlemail.com

Dependencies
------------

-  iconv
-  GNU patch


Special Update/Build Instructions
---------------------------------

-  TODO:
-  Perhaps also modify CXXFLAGS (and/or CPPFLAGS).
-  We currently don't use anything of GMP's/MPIR's CC setting, and
   matching with the current compiler (``$CC``) is perhaps suboptimal.
-  Remove some files / directories not needed for Sage from upstream:
-  build.vc\* directories (Microsoft Visual C build files)
-  3.0.0-644faf502c56f97d9accd301965fc57d6ec70868
   was created by running the spkg-src script.
