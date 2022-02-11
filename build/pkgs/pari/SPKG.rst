pari: Computer algebra system for fast computations in number theory
====================================================================

Description
-----------

PARI/GP is a widely used computer algebra system designed for fast
computations in number theory (factorizations, algebraic number theory,
elliptic curves...), but also contains a large number of other useful
functions to compute with mathematical entities such as matrices,
polynomials, power series, algebraic numbers etc., and a lot of
transcendental functions. PARI is also available as a C library to allow
for faster computations.

Originally developed by Henri Cohen and his co-workers (Université
Bordeaux I, France), PARI is now under the GPL and maintained by Karim
Belabas with the help of many volunteer contributors.

License
-------

GPL version 2+


Upstream Contact
----------------

-  http://pari.math.u-bordeaux.fr/

Dependencies
------------

-  Perl
-  MPIR or GMP
-  Readline
-  GNU patch (shipped with Sage)


Special Update/Build Instructions
---------------------------------

See patches/README.txt for a list of patches.

The current upstream tarball was created from the PARI git repository by
running "make snapshot".
