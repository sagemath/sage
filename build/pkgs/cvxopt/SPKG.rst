cvxopt: Python software for convex optimization
===============================================

Description
-----------

CVXOPT is a free software package for convex optimization based on the
Python programming language. It can be used with the interactive Python
interpreter, on the command line by executing Python scripts, or
integrated in other software via Python extension modules. Its main
purpose is to make the development of software for convex optimization
applications straightforward by building on Python's extensive standard
library and on the strengths of Python as a high-level programming
language.


Upstream Contact
----------------

-  J. Dahl <dahl.joachim@gmail.com>
-  L. Vandenberghe <vandenbe@ee.ucla.edu>

https://cvxopt.org/

License
-------

GPLv3 or later. Includes parts under GPLv2, GNU Lesser General Public
License, v2.1. See src/LICENSE for more details. (Sage-compatible)

Dependencies
------------

-  GNU patch
-  GSL
-  GLPK


Special Update/Build Instructions
---------------------------------

-  cvxopt.h.patch: Fix building with GCC on Solaris.

-  setup.py.patch: look for libraries and includes in $SAGE_LOCAL
   instead of /usr. Add fortran, blas,... libraries if needed.
   Build with GSL and GLPK support.

-  remove doc/html/, as it can be rebuild by invoking 'sage -sh' and
   running 'make html' in doc/

-  TODO: Add more tests in spkg-check

-  TODO: one might want to enhance the code to allow other Sage
   random sources, at the moment only GSL is used in CVXOPT-1.1.3
   spkg, apparently it will need an unclear to me "with seed(..)"
   construct.
