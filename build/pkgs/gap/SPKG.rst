gap: Groups, Algorithms, Programming - a system for computational discrete algebra
==================================================================================

Description
-----------

GAP is a system for computational discrete algebra, with particular
emphasis on Computational Group Theory. GAP provides a programming
language, a library of thousands of functions implementing algebraic
algorithms written in the GAP language as well as large data libraries
of algebraic objects. See also the overview and the description of the
mathematical capabilities. GAP is used in research and teaching for
studying groups and their representations, rings, vector spaces,
algebras, combinatorial structures, and more. The system, including
source, is distributed freely. You can study and easily modify or extend
it for your special use.

This is a stripped-down version of GAP. The databases, which are
architecture-independent, are in a separate package.


Upstream Contact
----------------

https://www.gap-system.org

Mailing list at https://mail.gap-system.org/mailman/listinfo/gap

Special Update/Build Instructions
---------------------------------

This is a stripped-down version of GAP. The downloading of the sources
and removal of unneeded parts is done by the script spkg-src. When you
update GAP, please also update and use the spkg-src script.

-  Do we really want to copy everything from the build directory???

   You need the full GAP tree to compile/install many GAP packages.

-  There's apparently a command missing (in ``spkg-install``) building
   the
   (HTML?) documentation. Earlier changelog entries as well as the
   description
   above state the documentation was removed from the upstream
   sources...
   Since the (pre-)built HTML documentation is currently included, I've
   commented out some lines in that part of ``spkg-install``. -leif

Patches
~~~~~~~
