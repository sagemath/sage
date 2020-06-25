Singular
========

Description
-----------

Singular is a computer algebra system for polynomial computations, with
special emphasis on commutative and non-commutative algebra, algebraic
geometry, and singularity theory.

License
-------

GPLv2 or GPLv3


Upstream Contact
----------------

libsingular-devel@mathematik.uni-kl.de

http://www.singular.uni-kl.de/

Dependencies
------------

-  GNU patch
-  readline
-  GMP/MPIR
-  MPFR
-  NTL
-  FLINT


Special Update/Build Instructions
---------------------------------

See spkg-src to create the source tarball.

Other notes:

-  If the environment variable SAGE_DEBUG is set to "yes", then
   omalloc will be replaced by xalloc. The resulting Singular executable
   and libsingular library will be slower than with omalloc, but allow
   for easier debugging of memory corruptions.
