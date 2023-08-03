p_group_cohomology: Modular cohomology rings of finite groups
=============================================================

Description
-----------

Modular Cohomology Rings of Finite Groups

The package is located at http://users.fmi.uni-jena.de/cohomology/,
that's to say the tarball p_group_cohomology-x.y.tar.xz can be found
there and the documentation of the package is provided at
http://users.fmi.uni-jena.de/cohomology/documentation/

License
-------

Copyright (C) 2018 Simon A. King <simon.king@uni-jena.de> Copyright (C)
2011 Simon A. King <simon.king@uni-jena.de> Copyright (C) 2009 Simon A.
King <simon.king@nuigalway.ie> and

   David J. Green <david.green@uni-jena.de>

Distributed under the terms of the GNU General Public License (GPL),
version 2 or later (at your choice).

   This code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.

The full text of the GPL is available at:

   http://www.gnu.org/licenses/

The package includes a data base of cohomology rings of the groups of
order 64 and provides access to a data base of cohomology rings of the
groups of order 128 and 243, located at

   http://cohomology.uni-jena.de/db/

These data bases are distributed under the Creative Commons
Attribution-Share Alike 3.0 License. The full text of this licence is
available at

   http://creativecommons.org/licenses/by-sa/3.0/


SPKG Maintainers
----------------

Simon A. King <simon.king@uni-jena.de>


Upstream Contact
----------------

Simon A. King <simon.king@uni-jena.de> David J. Green
<david.green@uni-jena.de>

Acknowledgements
----------------

The development of the initial version of this SPKG was funded by the
German Science Foundation, DFG project GR 1585/4.1, and was accomplished
at the Friedrich Schiller University Jena.

Since version 1.0.1, the further work on this SPKG was funded by Marie
Curie grant MTKD-CT-2006-042685 and was pursued at the National
University of Ireland, Galway. Since Novermber 2010, it is moved back to
Jena.

We thank William Stein for giving us access to various computers on
which we could build test the SPKG and on which some huge computations
could be completed, and acknowledge the support by National Science
Foundation Grant No. DMS-0821725.

We thank Mathieu Dutour Sikirić for hints on how to use GAP more
efficiently.

We owe Peter Symonds the idea of using the Poincaré series in a rather
efficient completeness criterion.

We are greatful to John Palmieri for his help on making
p_group_cohomology work with python-3.

Dependencies
------------

-  The SharedMeatAxe needs to be installed, as a build time dependency.

   This can be met by installing the meataxe spkg

Testing
-------

Our package provides a very short test suite for David Green's routines
for the computation of minimal projective resolutions. The majority of
this package's tests is formed by doc tests in the Cython code. In fact,
any class, method and function is covered by tests.

Note that internet access is required for these tests, as it is
attempted to download cohomology rings from a public data base in the
web.

The script ``spkg-check`` calls ``sage -t --force_lib`` on the files
in ``pGroupCohomology``.

Documentation
-------------

The documentation of this package is automatically built, if the
environment variable SAGE_SPKG_INSTALL_DOCS is yes (do "export
SAGE_SPKG_INSTALL_DOCS=yes" on the command line before installation).
The documents are put into
SAGE_ROOT/local/share/doc/p_group_cohomology/.
