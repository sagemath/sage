tides: Integration of ODEs
==========================

Description
-----------

TIDES is a library for integration of ODEs with high precision.

License
-------

GPLv3+


Upstream Contact
----------------

-  Marcos Rodriguez (marcos@unizar.es)

Dependencies
------------

-  gcc
-  mpfr
-  gmp


Special Update/Build Instructions
---------------------------------

minc_tides.patch changes the size of the name of the temporal files, so
there is no problem in systems that use long names. Also solves a bug in
the inverse function.
