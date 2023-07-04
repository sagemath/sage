.. -*- coding: utf-8 -*-

.. linkall

.. _tips:

=========================
Polyhedra tips and tricks
=========================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>


Operation shortcuts
=================================================

You can obtain different operations using natural symbols:

::

    sage: Cube = polytopes.cube()
    sage: Octahedron = 3/2*Cube.polar()  # Dilation
    sage: Cube + Octahedron   # Minkowski sum
    A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices
    sage: Cube & Octahedron   # Intersection
    A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices
    sage: Cube * Octahedron   # Cartesian product
    A 6-dimensional polyhedron in QQ^6 defined as the convex hull of 48 vertices
    sage: Cube - Polyhedron(vertices=[[-1,0,0],[1,0,0]])  # Minkowski difference
    A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices

.. end of output

Sage input function
==============================================================

If you are working with a polyhedron that was difficult to construct
and you would like to get back the proper Sage input code to reproduce this
object, you can!

::

    sage: Cube = polytopes.cube()
    sage: TCube = Cube.truncation().dilation(1/2)
    sage: sage_input(TCube)
    Polyhedron(backend='ppl', base_ring=QQ, vertices=[(1/6, -1/2, -1/2),
    (1/2, -1/6, -1/2), (1/2, 1/6, -1/2), (1/2, 1/2, -1/6), (1/2, 1/2, 1/6),
    (1/2, 1/6, 1/2), (1/6, 1/2, 1/2), (1/2, -1/6, 1/2), (1/6, 1/2, -1/2),
    (1/6, -1/2, 1/2), (1/2, -1/2, 1/6), (1/2, -1/2, -1/6), (-1/2, 1/6, -1/2),
    (-1/2, -1/2, 1/6), (-1/2, 1/6, 1/2), (-1/2, 1/2, 1/6), (-1/6, 1/2, 1/2),
    (-1/2, 1/2, -1/6), (-1/6, 1/2, -1/2), (-1/2, -1/6, 1/2), (-1/6, -1/2, 1/2),
    (-1/2, -1/2, -1/6), (-1/6, -1/2, -1/2), (-1/2, -1/6, -1/2)])

.. end of output


:code:`Hrepresentation_str`
==============================================================

If you would like to visualize the `H`-representation nicely and even get
the latex presentation, there is a method for that!

::

    sage: Nice_repr = TCube.Hrepresentation_str()
    sage: print(Nice_repr)
    -6*x0 - 6*x1 - 6*x2 >= -7
    -6*x0 - 6*x1 + 6*x2 >= -7
    -6*x0 + 6*x1 - 6*x2 >= -7
    -6*x0 + 6*x1 + 6*x2 >= -7
                  -2*x0 >= -1
                  -2*x1 >= -1
                  -2*x2 >= -1
     6*x0 + 6*x1 + 6*x2 >= -7
                   2*x2 >= -1
                   2*x1 >= -1
                   2*x0 >= -1
     6*x0 - 6*x1 - 6*x2 >= -7
     6*x0 - 6*x1 + 6*x2 >= -7
     6*x0 + 6*x1 - 6*x2 >= -7

    sage: print(TCube.Hrepresentation_str(latex=True))
    \begin{array}{rcl}
    -6 x_{0} - 6 x_{1} - 6 x_{2} & \geq & -7 \\
    -6 x_{0} - 6 x_{1} + 6 x_{2} & \geq & -7 \\
    -6 x_{0} + 6 x_{1} - 6 x_{2} & \geq & -7 \\
    -6 x_{0} + 6 x_{1} + 6 x_{2} & \geq & -7 \\
                        -2 x_{0} & \geq & -1 \\
                        -2 x_{1} & \geq & -1 \\
                        -2 x_{2} & \geq & -1 \\
     6 x_{0} + 6 x_{1} + 6 x_{2} & \geq & -7 \\
                         2 x_{2} & \geq & -1 \\
                         2 x_{1} & \geq & -1 \\
                         2 x_{0} & \geq & -1 \\
     6 x_{0} - 6 x_{1} - 6 x_{2} & \geq & -7 \\
     6 x_{0} - 6 x_{1} + 6 x_{2} & \geq & -7 \\
     6 x_{0} + 6 x_{1} - 6 x_{2} & \geq & -7
    \end{array}

    sage: Latex_repr = LatexExpr(TCube.Hrepresentation_str(latex=True))
    sage: view(Latex_repr)  # not tested

.. end of output

The `style` parameter allows to change the way to print the `H`-relations:

::

    sage: P = polytopes.permutahedron(3)
    sage: print(P.Hrepresentation_str(style='<='))
    -x0 - x1 - x2 == -6
         -x0 - x1 <= -3
          x0 + x1 <=  5
              -x1 <= -1
               x0 <=  3
              -x0 <= -1
               x1 <=  3
    sage: print(P.Hrepresentation_str(style='positive'))
    x0 + x1 + x2 == 6
         x0 + x1 >= 3
               5 >= x0 + x1
              x1 >= 1
               3 >= x0
              x0 >= 1
               3 >= x1

.. end of output
