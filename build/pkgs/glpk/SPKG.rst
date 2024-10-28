glpk: GNU Linear Programming Kit
================================

Description
-----------

The GLPK (GNU Linear Programming Kit) package is intended for solving
large-scale linear programming (LP), mixed integer programming (MIP),
and other related problems. It is a set of routines written in ANSI C
and organized in the form of a callable library.

GLPK supports the GNU MathProg modelling language, which is a subset of
the AMPL language.

The GLPK package includes the following main components:

-  primal and dual simplex methods
-  primal-dual interior-point method
-  branch-and-cut method
-  translator for GNU MathProg
-  application program interface (API)
-  stand-alone LP/MIP solver

License
-------

The GLPK package is GPL version 3.


Upstream Contact
----------------

GLPK is currently being maintained by:

-  Andrew Makhorin (mao@gnu.org, mao@mai2.rcnet.ru)

http://www.gnu.org/software/glpk/#maintainer

Special Update/Build Instructions
---------------------------------

-  ``configure`` doesn't support specifying the location of the GMP
   library to use; only ``--with-gmp[=yes]`` or ``--with-gmp=no``
   are valid options. (So we \*have to\* add Sage's include and
   library directories to ``CPPFLAGS`` and ``LDFLAGS``, respectively.)

-  Do we need the ``--disable-static``? The stand-alone solver presumably
   runs faster when built with a static library; also other
   (stand-alone)
   programs using it would.
   (Instead, we should perhaps use ``--enable-static --enable-shared``
   to
   go safe.)
