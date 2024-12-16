.. linkall

**********
Interfaces
**********

A central facet of Sage is that it supports computation with objects in
many different computer algebra systems "under one roof" using a
common interface and clean programming language.

The console and interact methods of an interface do very different
things. For example, using GAP as an example:


#. ``gap.console()``: This opens the GAP console - it transfers
   control to GAP. Here Sage is serving as nothing more than a
   convenient program launcher, similar to the Linux bash shell.

#. ``gap.interact()``: This is a convenient way to interact with a
   running GAP instance that may be "full of" Sage objects. You can
   import Sage objects into this GAP session (even from the
   interactive interface), etc.


.. index: PARI; GP

GP/PARI
=======

PARI is a compact, very mature, highly optimized C program whose
primary focus is number theory. There are two very distinct
interfaces that you can use in Sage:


-  ``gp`` -- the "**G** o **P** ARI" interpreter, and

-  ``pari`` -- the PARI C library.


For example, the following are two ways of doing the same thing.
They look identical, but the output is actually different, and what
happens behind the scenes is drastically different.

::

    sage: gp('znprimroot(10007)')
    Mod(5, 10007)
    sage: pari('znprimroot(10007)')
    Mod(5, 10007)

In the first case, a separate copy of the GP interpreter is started
as a server, and the string ``'znprimroot(10007)'`` is sent to it,
evaluated by GP, and the result is assigned to a variable in GP
(which takes up space in the child GP processes memory that won't
be freed). Then the value of that variable is displayed. In the
second case, no separate program is started, and the string
``'znprimroot(10007)'`` is evaluated by a certain PARI C library
function. The result is stored in a piece of memory on the Python
heap, which is freed when the variable is no longer referenced. The
objects have different types:

::

    sage: type(gp('znprimroot(10007)'))
    <class 'sage.interfaces.gp.GpElement'>
    sage: type(pari('znprimroot(10007)'))
    <class 'cypari2.gen.Gen'>

So which should you use? It depends on what you're doing. The GP
interface can do absolutely anything you could do in the usual
GP/PARI command line program, since it is running that program. In
particular, you can load complicated PARI programs and run them. In
contrast, the PARI interface (via the C library) is much more
restrictive. First, not all member functions have been implemented.
Second, a lot of code, e.g., involving numerical integration, won't
work via the PARI interface. That said, the PARI interface can be
significantly faster and more robust than the GP one.

(If the GP interface runs out of memory evaluating a given input
line, it will silently and automatically double the stack size and
retry that input line.  Thus your computation won't crash if you didn't
correctly anticipate the amount of memory that would be needed.  This
is a nice trick the usual GP interpreter doesn't seem to provide.  Regarding
the PARI C library interface, it immediately copies each created
object off of the PARI stack, hence the stack never grows.  However,
each object must not exceed 100MB in size, or the stack will overflow
when the object is being created.  This extra copying does impose
a slight performance penalty.)

In summary, Sage uses the PARI C library to provide functionality
similar to that provided by the GP/PARI interpreter, except with
different sophisticated memory management and the Python
programming language.

First we create a PARI list from a Python list.

::

    sage: v = pari([1,2,3,4,5])
    sage: v
    [1, 2, 3, 4, 5]
    sage: type(v)
    <class 'cypari2.gen.Gen'>

Every PARI object is of type ``Gen``. The PARI type of the
underlying object can be obtained using the ``type`` member
function.

::

    sage: v.type()
    't_VEC'

In PARI, to create an elliptic curve we enter
``ellinit([1,2,3,4,5])``. Sage is similar, except that ``ellinit`` is a
method that can be called on any PARI object, e.g., our
``t_VEC`` `v`.

::

    sage: e = v.ellinit()
    sage: e.type()
    't_VEC'
    sage: pari(e)[:13]
    [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351]

Now that we have an elliptic curve object, we can compute some
things about it.

::

    sage: e.elltors()
    [1, [], []]
    sage: e.ellglobalred()
    [10351, [1, -1, 0, -1], 1, [11, 1; 941, 1], [[1, 5, 0, 1], [1, 5, 0, 1]]]
    sage: f = e.ellchangecurve([1,-1,0,-1])
    sage: f[:5]
    [1, -1, 0, 4, 3]

.. index: GAP

.. _section-gap:

GAP
===

Sage comes with GAP for computational discrete mathematics,
especially group theory.

Here's an example of GAP's ``IdGroup`` function.

::

    sage: G = gap('Group((1,2,3)(4,5), (3,4))')
    sage: G
    Group( [ (1,2,3)(4,5), (3,4) ] )
    sage: G.Center()
    Group( () )
    sage: G.IdGroup()
    [ 120, 34 ]
    sage: G.Order()
    120

We can do the same computation in Sage without explicitly invoking the
GAP interface as follows:

::

    sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
    sage: G.center()
    Subgroup generated by [()] of (Permutation Group with generators [(3,4), (1,2,3)(4,5)])
    sage: G.group_id()
    [120, 34]
    sage: n = G.order(); n
    120

For some GAP functionality, you should install an optional
Sage packages. This can be done with the command::

    sage -i gap_packages


Singular
========


Singular provides a massive and mature library for Gröbner bases,
multivariate polynomial gcds, bases of Riemann-Roch spaces of a
plane curve, and factorizations, among other things. We illustrate
multivariate polynomial factorization using the Sage interface to
Singular (do not type the ``....:``):

::

    sage: R1 = singular.ring(0, '(x,y)', 'dp')
    sage: R1
    polynomial ring, over a field, global ordering
    //   coefficients: QQ
    //   number of vars : 2
    //        block   1 : ordering dp
    //                  : names    x y
    //        block   2 : ordering C
    sage: f = singular('9*y^8 - 9*x^2*y^7 - 18*x^3*y^6 - 18*x^5*y^6 +'
    ....:     '9*x^6*y^4 + 18*x^7*y^5 + 36*x^8*y^4 + 9*x^10*y^4 - 18*x^11*y^2 -'
    ....:     '9*x^12*y^3 - 18*x^13*y^2 + 9*x^16')

Now that we have defined :math:`f`, we print it and factor.

::

    sage: f
    9*x^16-18*x^13*y^2-9*x^12*y^3+9*x^10*y^4-18*x^11*y^2+36*x^8*y^4+18*x^7*y^5-18*x^5*y^6+9*x^6*y^4-18*x^3*y^6-9*x^2*y^7+9*y^8
    sage: f.parent()
    Singular
    sage: F = f.factorize(); F
    [1]:
       _[1]=9
       _[2]=x^6-2*x^3*y^2-x^2*y^3+y^4
       _[3]=-x^5+y^2
    [2]:
       1,1,2
    sage: F[1][2]
    x^6-2*x^3*y^2-x^2*y^3+y^4

As with the GAP example in :ref:`section-gap`, we can compute the
above factorization without explicitly using the Singular interface
(however, behind the scenes Sage uses the Singular interface for the
actual computation). Do not type the ``....:``:

::

    sage: x, y = QQ['x, y'].gens()
    sage: f = (9*y^8 - 9*x^2*y^7 - 18*x^3*y^6 - 18*x^5*y^6 + 9*x^6*y^4
    ....:     + 18*x^7*y^5 + 36*x^8*y^4 + 9*x^10*y^4 - 18*x^11*y^2 - 9*x^12*y^3
    ....:     - 18*x^13*y^2 + 9*x^16)
    sage: factor(f)
    (9) * (-x^5 + y^2)^2 * (x^6 - 2*x^3*y^2 - x^2*y^3 + y^4)

.. _section-maxima:

Maxima
======

Maxima is included with Sage, as well as a Lisp implementation. The
gnuplot package (which Maxima uses by default for plotting) is
distributed as a Sage optional package. Among other things, Maxima
does symbolic manipulation. Maxima can integrate and differentiate
functions symbolically, solve 1st order ODEs, most linear 2nd order
ODEs, and has implemented the Laplace transform method for linear ODEs
of any degree. Maxima also knows about a wide range of special
functions, has plotting capabilities via gnuplot, and has methods to
solve and manipulate matrices (such as row reduction, eigenvalues and
eigenvectors), and polynomial equations.

We illustrate the Sage/Maxima interface by constructing the matrix
whose :math:`i,j` entry is :math:`i/j`, for
:math:`i,j=1,\ldots,4`.

::

    sage: f = maxima.eval('ij_entry[i,j] := i/j')
    sage: A = maxima('genmatrix(ij_entry,4,4)'); A
    matrix([1,1/2,1/3,1/4],[2,1,2/3,1/2],[3,3/2,1,3/4],[4,2,4/3,1])
    sage: A.determinant()
    0
    sage: A.echelon()
    matrix([1,1/2,1/3,1/4],[0,0,0,0],[0,0,0,0],[0,0,0,0])
    sage: A.eigenvalues()
    [[0,4],[3,1]]
    sage: A.eigenvectors().sage()
    [[[0, 4], [3, 1]], [[[1, 0, 0, -4], [0, 1, 0, -2], [0, 0, 1, -4/3]], [[1, 2, 3, 4]]]]

Here's another example:

::

    sage: A = maxima("matrix ([1, 0, 0], [1, -1, 0], [1, 3, -2])")
    sage: eigA = A.eigenvectors()
    sage: V = VectorSpace(QQ,3)
    sage: eigA
    [[[-2,-1,1],[1,1,1]],[[[0,0,1]],[[0,1,3]],[[1,1/2,5/6]]]]
    sage: v1 = V(sage_eval(repr(eigA[1][0][0]))); lambda1 = eigA[0][0][0]
    sage: v2 = V(sage_eval(repr(eigA[1][1][0]))); lambda2 = eigA[0][0][1]
    sage: v3 = V(sage_eval(repr(eigA[1][2][0]))); lambda3 = eigA[0][0][2]

    sage: M = MatrixSpace(QQ,3,3)
    sage: AA = M([[1,0,0],[1, - 1,0],[1,3, - 2]])
    sage: b1 = v1.base_ring()
    sage: AA*v1 == b1(lambda1)*v1
    True
    sage: b2 = v2.base_ring()
    sage: AA*v2 == b2(lambda2)*v2
    True
    sage: b3 = v3.base_ring()
    sage: AA*v3 == b3(lambda3)*v3
    True

Finally, we give an example of using Sage to plot using ``openmath``.
Many of these were modified from the Maxima reference manual.

A 2D plot of several functions (do not type the ``....:``)::

    sage: maxima.plot2d('[cos(7*x),cos(23*x)^4,sin(13*x)^3]','[x,0,1]',  # not tested
    ....:     '[plot_format,openmath]')

A "live" 3D plot which you can move with your mouse (do not type
the ``....:``)::

    sage: maxima.plot3d ("2^(-u^2 + v^2)", "[u, -3, 3]", "[v, -2, 2]",  # not tested
    ....:     '[plot_format, openmath]')
    sage: maxima.plot3d("atan(-x^2 + y^3/4)", "[x, -4, 4]", "[y, -4, 4]",  # not tested
    ....:     "[grid, 50, 50]",'[plot_format, openmath]')

The next plot is the famous Möbius strip (do not type the ``....:``)::

    sage: maxima.plot3d("[cos(x)*(3 + y*cos(x/2)), sin(x)*(3 + y*cos(x/2)), y*sin(x/2)]",  # not tested
    ....:     "[x, -4, 4]", "[y, -4, 4]", '[plot_format, openmath]')

The next plot is the famous Klein bottle (do not type the ``....:``)::

    sage: maxima("expr_1: 5*cos(x)*(cos(x/2)*cos(y) + sin(x/2)*sin(2*y)+ 3.0) - 10.0")
    5*cos(x)*(sin(x/2)*sin(2*y)+cos(x/2)*cos(y)+3.0)-10.0
    sage: maxima("expr_2: -5*sin(x)*(cos(x/2)*cos(y) + sin(x/2)*sin(2*y)+ 3.0)").sage()
    -5*(cos(1/2*x)*cos(y) + sin(1/2*x)*sin(2*y) + 3.0)*sin(x)
    sage: maxima("expr_3: 5*(-sin(x/2)*cos(y) + cos(x/2)*sin(2*y))")
    5*(cos(x/2)*sin(2*y)-sin(x/2)*cos(y))
    sage: maxima.plot3d ("[expr_1, expr_2, expr_3]", "[x, -%pi, %pi]",  # not tested
    ....:     "[y, -%pi, %pi]", "['grid, 40, 40]", '[plot_format, openmath]')
