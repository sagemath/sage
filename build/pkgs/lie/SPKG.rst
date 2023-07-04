lie: Library for the representation theory of complex semisimple Lie groups and algebras
========================================================================================

Description
-----------

LiE is the name of a software package that enables mathematicians and
physicists to perform computations of a Lie group theoretic nature. It
focuses on the representation theory of complex semisimple (reductive)
Lie groups and algebras, and on the structure of their Weyl groups and
root systems.

LiE does not compute directly with elements of the Lie groups and
algebras themselves; it rather computes with weights, roots, characters
and similar objects. Some specialities of LiE are: tensor product
decompositions, branching to subgroups, Weyl group orbits, reduced
elements in Weyl groups, distinguished coset representatives and much
more. These operations have been compiled into the program which results
in fast execution: typically one or two orders of magnitude faster than
similar programs written in a general purpose program.

The LiE programming language makes it possible to customise and extend
the package with more mathematical functions. A user manual is provided
containing many examples.

LiE establishes an interactive environment from which commands can be
given that involve basic programming primitives and powerful built-in
functions. These commands are read by an interpreter built into the
package and passed to the core of the system. This core consists of
programs representing some 100 mathematical functions. The interpreter
offers on-line facilities which explain operations and functions, and
which give background information about Lie group theoretical concepts
and about currently valid definitions and values.

(from http://www-math.univ-poitiers.fr/~maavl/LiE/description.html )

License
-------

GNU Lesser General Public License (LGPL), version unspecified


Upstream Contact
----------------

-  Marc van Leeuwen, http://www-math.univ-poitiers.fr/~maavl/

Dependencies
------------

-  readline
-  ncurses
-  bison (not included in this package or in Sage!)
