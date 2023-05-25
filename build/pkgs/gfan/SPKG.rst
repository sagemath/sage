gfan: Groebner fans and tropical varieties
==========================================

Description
-----------

Gfan is a software package for computing Groebner fans and tropical
varieties.

These are polyhedral fans associated to polynomial ideals.
The maximal cones of a Groebner fan are in bijection with the marked
reduced Groebner bases of its defining ideal. The software computes all
marked reduced Groebner bases of an ideal. Their union is a universal
Groebner basis. The tropical variety of a polynomial ideal is a certain
subcomplex of the Groebner fan. Gfan contains algorithms for computing
this complex for general ideals and specialized algorithms for tropical
curves, tropical hypersurfaces and tropical varieties of prime ideals.
In addition to the above core functions the package contains many tools
which are useful in the study of Groebner bases, initial ideals and
tropical geometry. The full list of commands can be found in Appendix B
of the manual. For ordinary Groebner basis computations Gfan is not
competitive in speed compared to programs such as CoCoA, Singular and
Macaulay2.

License
-------

-  GPL version 2 or version 3 (according to the gfan website)


Upstream Contact
----------------

Anders Nedergaard Jensen

https://users-math.au.dk/jensen/software/gfan/gfan.html


Special Update/Build Instructions
---------------------------------

Remove the doc, homepage, and examples subdirectories, which take up
most of the space.
