coxeter3: Library for Coxeter groups, Bruhat ordering, Kazhdan-Lusztig polynomials
==================================================================================

Description
-----------

This package wraps Fokko Ducloux's Coxeter 3 C++ library

Features:

-  General Coxeter groups, implemented through the combinatorics of
   reduced words;
-  Reduced expression and normal form computations;
-  Bruhat ordering;
-  Ordinary Kazhdan-Lusztig polynomials;
-  Kazhdan-Lusztig polynomials with unequal parameters;
-  Inverse Kazhdan-Lusztig polynomials;
-  Cells and W-graphs;

http://math.univ-lyon1.fr/~ducloux/coxeter/coxeter3/english/coxeter3_e.html

This is a patched version done by Mike Hansen 2009-2013 and some fixes
by Nicolas M. Thiéry and Jean-Pierre Flori.

License
-------

GPL


Upstream Contact
----------------

github: https://github.com/tscrim/coxeter

Alas, Fokko Ducloux passed away in 2006.

http://math.univ-lyon1.fr/~ducloux/du_Cloux.html

Special Update/Build Instructions
---------------------------------

The source package was created by running ::

    commit=8ac9c71723c8ca57a836d6381aed125261e44e9e
    git clone https://github.com/tscrim/coxeter.git
    cd coxeter
    git archive $commit | bzip2 --best >coxeter-$commit.tar.bz2
