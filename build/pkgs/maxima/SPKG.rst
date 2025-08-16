maxima: System for manipulating symbolic and numerical expressions
==================================================================

Description
-----------

Maxima is a system for the manipulation of symbolic and numerical
expressions, including differentiation, integration, Taylor series,
Laplace transforms, ordinary differential equations, systems of linear
equations, polynomials, and sets, lists, vectors, matrices, and tensors.
Maxima yields high precision numeric results by using exact fractions,
arbitrary precision integers, and variable precision floating point
numbers. Maxima can plot functions and data in two and three dimensions.

For more information, see the Maxima web site

http://maxima.sourceforge.net

License
-------

Maxima is distributed under the GNU General Public License, with some
export restrictions from the U.S. Department of Energy. See the file
COPYING.


Upstream Contact
----------------

-  The Maxima mailing list - see
   http://maxima.sourceforge.net/maximalist.html

Special Update/Build Instructions
---------------------------------

1. Go to http://sourceforge.net/projects/maxima/files/Maxima-source/
   and download the source tarball maxima-x.y.z.tar.gz; place it in
   the upstream/ directory.

2. Update package-version.txt and run 'sage --package fix-checksum'.

3. Make sure the patches still apply cleanly, and update them if
   necessary.

4. Test the resulting package.

All patch files in the patches/ directory are applied. Descriptions of
these patches are either in the patch files themselves or below.

-  infodir.patch: Correct the path to the Info directory. Introduced
   in Issue #11348 (maxima test fails when install tree is moved).
