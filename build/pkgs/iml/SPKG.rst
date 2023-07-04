iml: Integer Matrix Library
===========================

Description
-----------

IML is a free library of C source code which implements algorithms for
computing exact solutions to dense systems of linear equations over the
integers. IML is designed to be used with the ATLAS/BLAS library and GMP
bignum library.

Written in portable C, IML can be used on both 32-bit and 64-bit
machines. It can be called from C++.

Website: https://www.cs.uwaterloo.ca/~astorjoh/iml.html

License
-------

-  GPLv2+


Upstream Contact
----------------

-  Zhuliang Chen z4chen@uwaterloo.ca
-  Arne Storjohann astorjoh@uwaterloo.ca

Special Update/Build Instructions
---------------------------------

-  As of version 1.0.4, you need to repackage the upstream tarball
   using the spkg-src script because there was a bugfix version of 1.0.4
   reposted upstream without version number bump.

Patches
~~~~~~~

-  examples.patch: Modified some of the examples.
