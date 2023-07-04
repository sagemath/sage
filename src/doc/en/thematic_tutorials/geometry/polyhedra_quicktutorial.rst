.. -*- coding: utf-8 -*-

.. linkall

.. _polytutorial:

A Brief Introduction to Polytopes in Sage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. MODULEAUTHOR:: sarah-marie belcastro <smbelcas@toroidalsnark.net>

If you already know some convex geometry  *a la*  Grünbaum or
Brøndsted, then you may have itched to get your hands dirty with some
polytope calculations.  Here is a mini\-guide to doing just that.

Basics
""""""

First, let's define a polytope as the convex hull of a set of points,
i.e. given  `S` we compute  `P={\rm conv}(S)`:


::

    sage: P1 = Polyhedron(vertices = [[-5,2], [4,4], [3,0], [1,0], [2,-4], [-3,-1], [-5,-3]])
    sage: P1
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices

.. end of output

Notice that Sage tells you the dimension of the polytope as well as the
dimension of the ambient space.

Of course, you want to know what this object looks like:


::

    sage: P1.plot()
    Graphics object consisting of 6 graphics primitives

.. end of output

Even in only 2 dimensions, it's a pain to figure out what the supporting
hyperplanes are.  Luckily Sage will take care of that for us.


::

    sage: for q in P1.Hrepresentation():
    ....:    print(q)
    An inequality (-4, 1) x + 12 >= 0
    An inequality (1, 7) x + 26 >= 0
    An inequality (1, 0) x + 5 >= 0
    An inequality (2, -9) x + 28 >= 0

.. end of output

That notation is not immediately parseable, because seriously,
those do not look like equations of lines (or of halfspaces, which is
really what they are).

``(-4, 1) x + 12 >= 0`` really means  `(-4, 1)\cdot\vec{x} + 12 \geq 0`.

So... if you want to define a polytope via inequalities, you have to
translate each inequality into a vector.  For example,
`(-4, 1)\cdot\vec{x} + 12 \geq 0` becomes (12, \-4, 1).


::

    sage: altP1 = Polyhedron(ieqs=[(12, -4, 1), (26, 1, 7),(5,1,0), (28, 2, -9)])
    sage: altP1.plot()
    Graphics object consisting of 6 graphics primitives

.. end of output

Other information you might want to pull out of Sage about a polytope is the
vertex list, which can be done in two ways:


::

    sage: for q in P1.Vrepresentation():
    ....:    print(q)
    A vertex at (-5, -3)
    A vertex at (-5, 2)
    A vertex at (4, 4)
    A vertex at (2, -4)

.. end of output

::

    sage: P1.vertices()
    (A vertex at (-5, -3), A vertex at (-5, 2), A vertex at (4, 4), A vertex at (2, -4))

.. end of output

Polar duals
"""""""""""

Surely you want to compute the polar dual:


::

    sage: P1dual = P1.polar()
    sage: P1dual
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

.. end of output

Check it out\-\-\-we started with an integer\-lattice polytope and dualized
to a rational\-lattice polytope.  Let's look at that.




::

    sage: P1dual.plot()
    Graphics object consisting of 6 graphics primitives


.. end of output

::

    sage: P1.plot() + P1dual.plot()
    Graphics object consisting of 12 graphics primitives


.. end of output

Oh, yeah, unless the polytope is unit\-sphere\-sized, the dual will be a
very different size.  Let's rescale.


::

    sage: ((1/4)*P1).plot() + (4*P1dual).plot()
    Graphics object consisting of 12 graphics primitives

.. end of output

If you think that looks a little bit shady, you're correct.  Here is an
example that makes the issue a bit clearer.


::

    sage: P2 = Polyhedron(vertices = [[-5,0], [-1,1], [-2,0], [1,0], [-2,-1], [-3,-1], [-5,-1]])
    sage: P2
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 5 vertices
    sage: P2dual = P2.polar(); P2dual
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 5 vertices
    sage: P2.plot() + P2dual.plot()
    Graphics object consisting of 14 graphics primitives

.. end of output

That is clearly not computing what we think of as the polar dual.  But look
at this...


::

    sage: P2.plot() + (-1*P2dual).plot()
    Graphics object consisting of 14 graphics primitives

.. end of output

Here is what's going on.

If a polytope ``P`` is in `\ZZ`, then...

(1) ...the dual is inverted in some way, which is vertically for polygons.

(2) ...the dual is taken of P itself.

(3) ...if the origin is not in P, then an error is returned.

However, if a polytope is  *not*  in `\ZZ`, for example if it's in `\QQ` or
``RDF``, then...

(1') ...the dual is not inverted.

(2') ...the dual is taken of P\-translated\-so\-barycenter\-is\-at\-origin.

Keep all of this in mind as you take polar duals.



Polytope Constructions
""""""""""""""""""""""

Minkowski sums!  Now with two syntaxes!


::

    sage: P1+P2
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 8 vertices

.. end of output

::

    sage: P1.minkowski_sum(P2)
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 8 vertices

.. end of output

Okay, fine.  We should have some 3\-dimensional examples, at least.
(Note that in order to display polytopes effectively you'll need
visualization software such as Javaview and Jmol installed.)


::

    sage: P3 = Polyhedron(vertices=[(0,0,0), (0,0,1/2), (0,1/2,0), (1/2,0,0), (3/4,1/5,3/2)]); P3
    A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 5 vertices
    sage: P4 = Polyhedron(vertices=[(-1,1,0),(1,1,0),(-1,0,1), (1,0,1),(0,-1,1),(0,1,1)]); P4
    A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
    sage: P3.plot() + P4.plot()
    Graphics3d Object

.. end of output

::

    sage: (P3+P4).plot()
    Graphics3d Object

.. end of output

We can also find the intersection of two polytopes... and this too has two
syntaxes!


::

    sage: int12 = P1.intersection(P2*.5); int12.plot()
    Graphics object consisting of 7 graphics primitives

.. end of output

::

    sage: int34 = P3 & P4; int34.plot()
    Graphics3d Object

.. end of output

Should one wish to translate, one can.


::

    sage: transP2 = P2.translation([2,1])
    sage: P2.plot() + transP2.plot()
    Graphics object consisting of 14 graphics primitives

.. end of output

Then of course we can take prisms, pyramids, and bipyramids of polytopes...


::

    sage: P2.prism().plot()
    Graphics3d Object

.. end of output

::

    sage: P1.pyramid().plot()
    Graphics3d Object

.. end of output

::

    sage: P2dual.bipyramid().plot()
    Graphics3d Object

.. end of output

Okay, fine.  Yes, Sage has some kinds of polytopes built in.
If you type ``polytopes.`` and then press ``TAB`` after the period, you'll get a
list of pre\-built polytopes.


::

    sage: P5 = polytopes.hypercube(5)
    sage: P6 = polytopes.cross_polytope(3)
    sage: P7 = polytopes.simplex(7)


.. end of output

Let's look at a 4\-dimensional polytope.


::

    sage: P8 = polytopes.hypercube(4)
    sage: P8.plot()
    Graphics3d Object

.. end of output

We can see it from a different perspective:


::

    sage: P8.schlegel_projection(position=1/2).plot()
    Graphics3d Object

.. end of output

Queries to polytopes
""""""""""""""""""""

Once you've constructed some polytope, you can ask Sage questions about it.


::

    sage: P1.contains([1,0])
    True

.. end of output

::

    sage: P1.interior_contains([3,0])
    False

.. end of output

::

    sage: P3.contains([1,0,0])
    False

.. end of output

Face information can be useful.


::

    sage: int34.f_vector()
    (1, 8, 12, 6, 1)

.. end of output

Well, geometric information might be  *more*  helpful...
Here we are told which of the vertices form each 2\-face:


::

    sage: [f.ambient_V_indices() for f in int34.faces(2)]
    [(2, 6, 7), (0, 1, 3, 5), (1, 3, 4), (0, 5, 6, 7), (0, 1, 2, 4, 6), (2, 3, 4, 5, 7)]

.. end of output

Yeah, that isn't so useful as it is.  Let's figure out the vertex and
hyperplane representations of the first face in the list.


::

    sage: first2faceofint34 = int34.faces(2)[0]
    sage: first2faceofint34.ambient_Hrepresentation(); first2faceofint34.vertices()
    (An inequality (0, 0, -1) x + 1 >= 0,)
    (A vertex at (2/3, 2/15, 1), A vertex at (3/8, 1/10, 1), A vertex at (1/2, 3/10, 1))

.. end of output

If you want more... :ref:`sage.geometry.polyhedron.base` is the first place you want to go.
