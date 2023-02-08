sympow: Computes special values of symmetric power elliptic curve L-functions
=============================================================================

Description
-----------

SYMPOW is a package to compute special values of symmetric power
elliptic curve L-functions. It can compute up to about 64 digits of
precision.

License
-------

-  See the file src/COPYING


Upstream Contact
----------------

SYMPOW does not appear to be maintained any longer.
Mark Watkins, the package author, now works at Magma.
Previous (possibly still usable) email is watkins@maths.usyd.edu.au

New upstream: https://gitlab.com/rezozer/forks/sympow

Dependencies
------------

-  GNU patch


Special Update/Build Instructions
---------------------------------

-  Some of the code is very dubious, and it is anyones guess really what
   the compiler does with it. For example, the following line exists in
   src/eulerfactors.c:

   if ((HECKE) && (d==1)) return hecke_good(p,ap,m,v);

   But since hecke_good is defined as returning void, it's hard to know
   exactly how this code behaves. I would not be surprised by any bugs
   that might show up. I (David Kirkby) would personally not trust this
   code much at all.

-  This is a difficult package to maintain. A github issue (#9758) has
   been
   opened to implement Watkins-Delaunay's algorithm for computing
   modular
   degrees in Sage. Once implemented, it should be possible to remove
   this
   package.

-  The package is configured such that the data files are in a directory
   below where 'sympow' is installed. If Sage is installed globally,
   then
   it will be impossible to create the data files without being root.
   This has been fixed in the Gentoo Linux distribution. Some
   information
   from Christopher can be seen on
   https://github.com/sagemath/sage/issues/9703
   This package will generate binary versions of all shipped datafiles,
   so these will work. However, creating totally new datafiles from
   scratch
   will not work.
