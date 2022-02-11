scipoptsuite: Mixed integer programming solver
==============================================

Description
-----------

SCIP is currently one of the fastest non-commercial mixed integer
programming (MIP) solvers. It is also a framework for constraint integer
programming and branch-cut-and-price. It allows total control of the
solution process and the access of detailed information down to the guts
of the solver.

License
-------

ZIB Academic License

The ZIB Academic License allows the use of software distributed under
this license without charge for research purposes as a member of a
non-commercial and academic institution, e.g., a university. The
software is available with its source code.

http://scip.zib.de/academic.txt


SPKG Maintainers
----------------

-  Martin Albrecht (original spkg)
-  Matthias Koeppe (updates for new spkg style)


Upstream Contact
----------------

http://scip.zib.de/doc/html/AUTHORS.shtml

Dependencies
------------

cmake


Special Update/Build Instructions
---------------------------------

We do not have permission to redistribute SCIP or SoPlex. Hence, you
must download it yourself from http://scip.zib.de and put the tarball
``scipoptsuite-VERSION.tgz`` in ``$SAGE_ROOT/upstream``, renaming
it to ``scipoptsuite-VERSION-do-not-distribute.tgz``.
