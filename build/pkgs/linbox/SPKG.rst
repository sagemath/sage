LinBox
======

Description
-----------

From http://linalg.org/: LinBox is a C++ template library for exact,
high-performance linear algebra computation with dense, sparse, and
structured matrices over the integers and over finite fields.

License
-------

LGPL V2 or later


Upstream Contact
----------------

-  <linbox-devel@googlegroups.com>
-  <linbox-use@googlegroups.com>


SPKG Repository
---------------

   https://bitbucket.org/malb/linbox-spkg

Dependencies
------------

-  GNU patch
-  GMP/MPIR
-  MPFR
-  NTL
-  fpLLL
-  IML
-  M4RI
-  M4RIE
-  Givaro
-  FFLAS/FFPACK
-  ATLAS (non-OSX)/The Accelerate FrameWork (on OSX)
-  ATLAS (non-MacOS X) / The Accelerate FrameWork (on MacOS X), or GSL's
   CBLAS


Special Update/Build Instructions
---------------------------------

TODO:

-  spkg-check is disabled for now, should work in the next release
   after 1.3.2.

-  Check whether ``make fullcheck`` works/builds, is worth running, and
   doesn't
   take ages. (Version 1.1.6 doesn't seem to have such a target.)
