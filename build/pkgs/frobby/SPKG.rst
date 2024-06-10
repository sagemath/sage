frobby: Computations on monomial ideals
=======================================

Description
-----------

The software package Frobby provides a number of computations on
monomial ideals. The current main feature is the socle of a monomial
ideal, which is largely equivalent to computing the maximal standard
monomials, the Alexander dual or the irreducible decomposition.

Operations on monomial ideals are much faster than algorithms designed
for ideals in general, which is what makes a specialized library for
these operations on monomial ideals useful.

License
-------

-  GPL version 2.0 or later


Upstream Contact
----------------

- http://www.broune.com/frobby/

- https://github.com/Macaulay2/frobby

Special Update/Build instructions
---------------------------------

Download Frobby at www.broune.com/ and then type "make spkg VER=blah"
which wil create an spkg named frobby-VER.spkg in bin/. The files
related to doing this is in the sage/ sub-directory of the Frobby source
distribution.
