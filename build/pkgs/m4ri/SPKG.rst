m4ri: fast arithmetic with dense matrices over GF(2)
====================================================

Description
-----------

M4RI: Library for matrix multiplication, reduction and inversion over
GF(2). (See also m4ri/README for a brief overview.)

License
-------

-  GNU General Public License Version 2 or later (see src/COPYING)


Upstream Contact
----------------

-  Authors: Martin Albrecht et al.
-  Email: <m4ri-devel@googlegroups.com>
-  Website: https://bitbucket.org/malb/m4ri

Special Update/Build Instructions
---------------------------------

-  Delete the upstream Mercurial repositories (file m4ri/.hgtags,
   directory m4ri/.hg).
-  Delete the directory m4ri/autom4te.cache (if present).
-  Delete m4ri.vcproj (and perhaps other unnecessary baggage).
-  Touch m4ri/configure to make sure it is newer than its sources.
