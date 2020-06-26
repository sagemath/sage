symmetrica
==========

Description
-----------

Symmetrica is a program developed by Lehrstuhl Mathematik II of the
University of Bayreuth. It has routines to handle the following topics

-  ordinary representation theory of the symmetric group and related
   groups (2/11/04)
-  ordinary representation theory of the classical groups
-  modular representation theory of the symmetric group
-  projective representation theory of the symmetric group
-  combinatorics of tableaux
-  symmetric functions and polynomials (7/22/04)
-  commutative and non commutative Schubert polynomials
-  operations of finite groups.
-  ordinary representation theory of Hecke algebras of type A_n

For more details check http://www.symmetrica.de (currently redirects to
http://www.algorithm.uni-bayreuth.de/en/research/SYMMETRICA)

License
-------

Public Domain (see the above web site)


Upstream Contact
----------------

-  Axel Kohnert - see http://www.mathe2.uni-bayreuth.de/axel/

Dependencies
------------

-  GNU patch (for applying the patches to upstream)


Special Update/Build Instructions
---------------------------------

The following patches are applied in spkg-install:

-  ``bruch.patch``: store integers in a temporary variable before freeing
   memory
-  ``de.patch``: turn off banner
-  ``int32.patch``: use ``int32_t`` and ``uint32_t`` for type INT.
-  ``sort_sum_rename.patch``: rename ``sort`` to ``sym_sort``, ``sum`` to ``sym_sum``
-  We copy over our own ``Makefile``:

   ``patches/makefile`` (Fix compiler, i.e., use ``$CC``, and let it use
   ``$CFLAGS``.)

Permissions in the upstream tarball are funky, please run ``chmod 644
src/*`` after unpacking.
