singular: Computer algebra system for polynomial computations, algebraic geometry, singularity theory
=====================================================================================================

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

https://www.singular.uni-kl.de/

Special Update/Build Instructions
---------------------------------

Other notes:

-  If the environment variable SAGE_DEBUG is set to "yes", then
   omalloc will be replaced by xalloc. The resulting Singular executable
   and libsingular library will be slower than with omalloc, but allow
   for easier debugging of memory corruptions.
