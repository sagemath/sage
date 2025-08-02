semigroups: An optional GAP package
===================================

Description
-----------

Installing this SPKG will install the corresponding GAP package, but
before you can use them in Sage, they still have to be loaded into
either the GAP interface or libgap::

  sage: gap.eval('LoadPackage("semigroups")')  # optional - semigroups
  'true'
  sage: libgap.LoadPackage("semigroups")       # optional - semigroups
  true

Those correspond to::

  gap> LoadPackage("semigroups");

within the GAP interface and libgap, respectively.

Upstream Contact
----------------

See https://semigroups.github.io/Semigroups/

Dependencies
------------

-  GAP (a standard spkg), gap_packages and libsemigroups (optional packages)

Notes
-----------
This is a GAP package for semigroups, and monoids. There are
particularly efficient methods for finitely presented semigroups and monoids,
and for semigroups and monoids consisting of transformations, partial
permutations, bipartitions, partitioned binary relations, subsemigroups of
regular Rees 0-matrix semigroups, and matrices of various semirings including
boolean matrices, matrices over finite fields, and certain tropical matrices.
Semigroups contains efficient methods for creating semigroups, monoids, and
inverse semigroups and monoids, calculating their Green's structure, ideals,
size, elements, group of units, small generating sets, testing membership,
finding the inverses of a regular element, factorizing elements over the
generators, and so on. It is possible to test if a semigroup satisfies a
particular property, such as if it is regular, simple, inverse, completely
regular, and a large number of further properties. There are methods for
finding presentations for a semigroup, the congruences of a semigroup, the
maximal subsemigroups of a finite semigroup, smaller degree partial
permutation representations, and the character tables of inverse semigroups.
There are functions for producing pictures of the Green's structure of a
semigroup, and for drawing graphical representations of certain types of
elements.
(Authors: James Mitchell, Marina Anagnostopoulou-Merkouri,
Thomas Breuer, Stuart Burrell, Reinis Cirpons, Tom Conti-Leslie,
Joseph Edwards, Attila Egri-Nagy, Luke Elliott, Fernando Flores Brito,
Tillman Froehlich, Nick Ham, Robert Hancock, Max Horn, Christopher Jefferson,
Julius Jonusas, Chinmaya Nagpal, Olexandr Konovalov, Artemis Konstantinidi,
Hyeokjun Kwon, Dima V. Pasechnik, Markus Pfeiffer, Christopher Russell,
Jack Schmidt, Sergio Siccha, Finn Smith, Ben Spiers, Nicolas Thiéry,
Maria Tsalakou, Chris Wensley, Murray Whyte, Wilf A. Wilson, Tianrun Yang,
Michael Young and Fabian Zickgraf)
