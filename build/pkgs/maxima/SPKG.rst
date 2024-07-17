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

-  0001-taylor2-Avoid-blowing-the-stack-when-diff-expand-isn.patch:
   Fix for Maxima bug #2520 (abs_integrate fails on abs(sin(x)) and
   abs(cos(x))). Introduced in Issue #13364 (Upgrade Maxima to
   5.29.1).

-  build-fasl.patch: Build a fasl library for ecl in addition to an
   executable program. Introduced in Issue #16178 (Build maxima fasl
   without asdf).

-  infodir.patch: Correct the path to the Info directory. Introduced
   in Issue #11348 (maxima test fails when install tree is moved).

-  matrixexp.patch: Fix matrixexp(matrix([%i*%pi])), which broke after
   Maxima 5.29.1. Introduced in Issue #13973.

-  maxima.system.patch: Set ``c::*compile-in-constants*`` to t.
   Introduced in Issue #11966 (OS X 10.7 Lion: Maxima fails to build).

-  undoing_true_false_printing_patch.patch: Revert an upstream change
   causing '?' to be printed around some words. Introduced in Trac
   #13364 (Upgrade Maxima to 5.29.1).
