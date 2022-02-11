.. -*- coding: utf-8 -*-

.. linkall

.. _polyhedron_visualization:

==================================================
Visualization of polyhedron objects in Sage
==================================================

.. MODULEAUTHOR:: sarah-marie belcastro <smbelcas@toroidalsnark.net>, Jean-Philippe Labbé <labbe@math.fu-berlin.de>


There are different ways to visualize polyhedron object of dimension at most 4.

:code:`render_solid`
==================================================

This plots the polyhedron as a solid. You can also adjust the :code:`opacity`
parameter.

::

    sage: Cube = polytopes.cube()
    sage: Cube.render_solid(opacity=0.7)
    Graphics3d Object

.. end of output

:code:`render_wireframe`
==================================================

This plots the graph (with unbounded edges) of the polyhedron

::

    sage: Cube.render_wireframe()
    Graphics3d Object

.. end of output

:code:`plot`
==================================================

The :code:`plot` method draws the graph, the polygons and vertices of the
polyhedron all together.

::

    sage: Cube.plot()
    Graphics3d Object

.. end of output

:code:`show`
==================================================

This is similar to :code:`plot` but does not return an object that you can
manipulate.


:code:`schlegel_projection`
==================================================

It is possible to visualize 4-dimensional polytopes using a schlegel diagram.

::

    sage: HC = polytopes.hypercube(4)
    sage: HC.schlegel_projection()
    The projection of a polyhedron into 3 dimensions
    sage: HC.schlegel_projection().plot()
    Graphics3d Object

.. end of output

We can see it from a different perspective by choosing point at a different
distance:

::

    sage: HC.schlegel_projection(position=1/4).plot()
    Graphics3d Object

.. end of output

It is possible to choose from which facet one sees the projection:

::

    sage: tHC = HC.face_truncation(HC.faces(0)[0])
    sage: tHC.facets()
    (A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 10 vertices,
     A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 10 vertices,
     A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 10 vertices,
     A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 10 vertices,
     A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 4 vertices,
     A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 8 vertices,
     A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 8 vertices,
     A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 8 vertices,
     A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 8 vertices)
    sage: tHC.schlegel_projection(tHC.facets()[4]).plot()
    Graphics3d Object

.. end of output

:code:`tikz`
==================================================

This method returns a tikz picture of the polytope (must be 2 or
3-dimensional). For more detail see the tutorial :ref:`polytikz`.

::

    sage: c = polytopes.cube()
    sage: c.tikz().splitlines()[:5]
    ['\\begin{tikzpicture}%',
     '\t[x={(1.000000cm, 0.000000cm)},',
     '\ty={(-0.000000cm, 1.000000cm)},',
     '\tz={(0.000000cm, -0.000000cm)},',
     '\tscale=1.000000,']

.. end of output
