******************
Algebraic Geometry
******************

.. index::
   pair: elliptic curve; point counting

Point counting on curves
========================
How do you count points on an elliptic curve over a finite field in
Sage?

Over prime finite fields, includes both the baby step giant step
method and the SEA (Schoof-Elkies-Atkin) algorithm (implemented in PARI
by Christophe Doche and Sylvain Duquesne). An example taken form the
Reference manual:

::

    sage: E = EllipticCurve(GF(10007),[1,2,3,4,5])
    sage: E.cardinality()
    10076

The command ``E.points()`` will return the actual list of rational
points.

How do you count points on a plane curve over a finite field? The
``rational_points`` command produces points by a simple enumeration
algorithm. Here is an example of the syntax:

::

    sage: x,y,z = PolynomialRing(GF(5), 3, 'xyz').gens()
    sage: C = Curve(y^2*z^7 - x^9 - x*z^8); C
    Projective Plane Curve over Finite Field of size 5 defined by -x^9 + y^2*z^7 - x*z^8
    sage: C.rational_points()
    [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1), (3 : 1 : 1), (3 : 4 : 1)]
    sage: C.rational_points(algorithm="bn")
    [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1), (3 : 1 : 1), (3 : 4 : 1)]

The option ``algorithm="bn`` uses Sage's Singular interface and
calls the ``brnoeth`` package.

Here is another example using Sage's ``rational_points`` applied to
Klein's quartic over :math:`GF(8)`.

::

    sage: x, y, z = PolynomialRing(GF(8,'a'), 3, 'xyz').gens()
    sage: f = x^3*y+y^3*z+x*z^3
    sage: C = Curve(f); C
    Projective Plane Curve over Finite Field in a of size 2^3 defined by x^3*y + y^3*z + x*z^3
    sage: C.rational_points()
    [(0 : 0 : 1),
     (0 : 1 : 0),
     (1 : 0 : 0),
     (1 : a : 1),
     (1 : a^2 : 1),
     (1 : a^2 + a : 1),
     (a : 1 : 1),
     (a : a^2 : 1),
     (a : a^2 + 1 : 1),
     (a + 1 : a + 1 : 1),
     (a + 1 : a^2 : 1),
     (a + 1 : a^2 + a + 1 : 1),
     (a^2 : 1 : 1),
     (a^2 : a^2 + a : 1),
     (a^2 : a^2 + a + 1 : 1),
     (a^2 + 1 : a + 1 : 1),
     (a^2 + 1 : a^2 + 1 : 1),
     (a^2 + 1 : a^2 + a : 1),
     (a^2 + a : 1 : 1),
     (a^2 + a : a : 1),
     (a^2 + a : a + 1 : 1),
     (a^2 + a + 1 : a : 1),
     (a^2 + a + 1 : a^2 + 1 : 1),
     (a^2 + a + 1 : a^2 + a + 1 : 1)]

Other methods
-------------


-  For a plane curve, you can use Singular's ``closed_points``
   command. The input is the vanishing ideal :math:`I` of the curve
   :math:`X` in a ring of :math:`2` variables :math:`F[x,y]`.
   The ``closed_points`` command returns a list of prime ideals (each a
   Gröbner basis), corresponding to the (distinct affine closed)
   points of :math:`V(I)`. Here's an example:

   .. skip

   ::

       sage: singular_console()
                            SINGULAR                             /  Development
        A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                              0<
            by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
       FB Mathematik der Universitaet, D-67653 Kaiserslautern    \
       // ** executing /home/wdj/sagefiles/sage-0.9.4/local/LIB/.singularrc
       > LIB "brnoeth.lib";
       > ring s = 2,(x,y),lp;
       > ideal I = x4+x,y4+y;
       > list L = closed_points(I);
       > L;
       [1]:
          _[1] = y
          _[2] = x
       [2]:
          _[1] = y
          _[2] = x+1
       [3]:
          _[1] = y
          _[2] = x2+x+1
       [4]:
          _[1] = y+1
          _[2] = x
       [5]:
          _[1] = y+1
          _[2] = x+1
       [6]:
          _[1] = y+1
          _[2] = x2+x+1
       [7]:
          _[1] = y2+y+1
          _[2] = x+1
       [8]:
          _[1] = y2+y+1
          _[2] = x
       [9]:
          _[1] = y2+y+1
          _[2] = x+y
       [10]:
          _[1] = y2+y+1
          _[2] = x+y+1
       > Auf Wiedersehen.

   ::

       sage: singular.lib("brnoeth.lib")
       sage: s = singular.ring(2,'(x,y)','lp')
       sage: I = singular.ideal('x^4+x', 'y^4+y')
       sage: L = singular.closed_points(I)
       sage: # Here you have all the points :
       sage: L       # random
       [1]:
          _[1]=y+1
          _[2]=x+1
       ...
       sage: l=[L[k].sage() for k in [1..10]]; len(l) # there are 10 points
       10
       sage: r=sorted(l[0].ring().gens()); r
       [y, x]
       sage: r in [t.gens() for t in l] #  one of them is given by [y,x]
       True

-  Another way to compute rational points is to use Singular's
   ``NSplaces`` command. Here's the Klein quartic over :math:`GF(8)`
   done this way:

   ::

       sage: singular.LIB("brnoeth.lib")
       sage: s = singular.ring(2,'(x,y)','lp')
       sage: f = singular.poly('x3y+y3+x')
       sage: klein1 = f.Adj_div(); print(klein1)
       [1]:
          [1]:
             // coefficients: ZZ/2...
       // number of vars : 2
       //        block   1 : ordering lp
       //                  : names    x y
       //        block   2 : ordering C
       ...
       sage: # define a curve X = {f = 0} over GF(2)
       sage: klein2 = singular.NSplaces(3,klein1)
       sage: print(singular.eval('extcurve(3,%s)'%klein2.name()))
       Total number of rational places : NrRatPl = 23
       ...
       sage: klein3 = singular.extcurve(3, klein2)

   Above we defined a curve :math:`X = \{f = 0\}` over
   :math:`GF(8)` in Singular.

   .. link

   ::

       sage: print(klein1)
       [1]:
          [1]:
             // coefficients: ZZ/2...
       // number of vars : 2
       //        block   1 : ordering lp
       //                  : names    x y
       //        block   2 : ordering C
          [2]:
             // coefficients: ZZ/2...
       // number of vars : 3
       //        block   1 : ordering lp
       //                  : names    x y z
       //        block   2 : ordering C
       [2]:
          4,3
       [3]:
          [1]:
             1,1
          [2]:
             1,2
       [4]:
          0
       [5]:
          [1]:
             [1]:
                // coefficients: ZZ/2...
       // number of vars : 3
       //        block   1 : ordering ls
       //                  : names    x y t
       //        block   2 : ordering C
             [2]:
                1,1
       sage: print(klein1[3])
       [1]:
          1,1
       [2]:
          1,2

   For the places of degree :math:`3`:

   .. link

   ::

       sage: print(klein2[3])
       [1]:
          1,1
       [2]:
          1,2
       [3]:
          3,1
       [4]:
          3,2
       [5]:
          3,3
       [6]:
          3,4
       [7]:
          3,5
       [8]:
          3,6
       [9]:
          3,7

   Each point below is a pair: (degree, point index number).

   .. link

   ::

       sage: print(klein3[3])
       [1]:
          1,1
       [2]:
          1,2
       [3]:
          3,1
       [4]:
          3,2
       [5]:
          3,3
       [6]:
          3,4
       [7]:
          3,5
       [8]:
          3,6
       [9]:
          3,7

   To actually get the points of :math:`X(GF(8))`:

   .. link

   ::

       sage: R = klein3[1][5]
       sage: R.set_ring()
       sage: singular("POINTS;")
       [1]:
          [1]:
             0
          [2]:
             1
          [3]:
             0
       [2]:
          [1]:
             1
          [2]:
             0
          [3]:
             0
       ...

   plus 21 others (omitted). There are a total of :math:`23`
   rational points.

.. index:: Riemann-Roch space

Riemann-Roch spaces using Singular
==================================

To compute a basis of the Riemann-Roch space of a divisor :math:`D`
on a curve over a field :math:`F`, one can use Sage's wrapper
``riemann_roch_basis`` of Singular's implementation of the Brill
Noether algorithm. Note that this wrapper currently only works when
:math:`F` is prime and the divisor :math:`D` is supported on rational points.
Below are examples of how to use ``riemann_roch_basis`` and how to use
Singular itself to help an understanding of how the wrapper works.

-  Using ``riemann_roch_basis``:

   ::

       sage: x, y, z = PolynomialRing(GF(5), 3, 'xyz').gens()
       sage: f = x^7 + y^7 + z^7
       sage: X = Curve(f); pts = X.rational_points()
       sage: D = X.divisor([ (3, pts[0]), (-1,pts[1]), (10, pts[5]) ])
       sage: X.riemann_roch_basis(D)
       [(-2*x + y)/(x + y), (-x + z)/(x + y)]

-  Using Singular's ``BrillNoether`` command (for details see the section
   Brill-Noether in the Singular online documentation
   (https://www.singular.uni-kl.de/Manual/4-3-0/sing_2254.htm and the
   paper {CF}):

   ::

       sage: singular.LIB('brnoeth.lib')
       sage: _ = singular.ring(5,'(x,y)','lp')
       sage: print(singular.eval("list X = Adj_div(-x5+y2+x);"))
       Computing affine singular points ...
       Computing all points at infinity ...
       Computing affine singular places ...
       Computing singular places at infinity ...
       Computing non-singular places at infinity ...
       Adjunction divisor computed successfully
       <BLANKLINE>
       The genus of the curve is 2
       sage: print(singular.eval("X = NSplaces(1,X);"))
       Computing non-singular affine places of degree 1 ...
       sage: print(singular("X[3];"))
       [1]:
          1,1
       [2]:
          1,2
       [3]:
          1,3
       [4]:
          1,4
       [5]:
          1,5
       [6]:
          1,6

   The first integer of each pair in the above list is the degree
   `d` of a point. The second integer is the index of this point
   in the list POINTS of the ring X[5][`d`][1]. Note that the
   order of this latter list is different every time the algorithm
   is run, e.g. `1`, `1` in the above list refers to a different
   rational point each time. A divisor is given by defining a list
   `G` of integers of the same length as X[3] such that if the
   `k`-th entry of X[3] is `d`, `i`, then the `k`-th entry of `G` is
   the multiplicity of the divisor at the `i`-th point in the list
   POINTS of the ring X[5][`d`][1]. Let us proceed by defining a
   "random" divisor of degree 12 and computing a basis of its
   Riemann-Roch space:

   .. link

   ::

       sage: singular.eval("intvec G = 4,4,4,0,0,0;")
       ''
       sage: singular.eval("def R = X[1][2];")
       ''
       sage: singular.eval("setring R;")
       ''
       sage: print(singular.eval("list LG = BrillNoether(G,X);"))
       Forms of degree 6 :
       28
       <BLANKLINE>
       Vector basis successfully computed
       <BLANKLINE>


.. index::
   pair: codes; algebraic-geometric

AG codes
--------

Sage can compute an AG code :math:`C=C_X(D,E)` by calling
Singular's BrillNoether to compute a basis of the Riemann Roch
space :math:`L(D)=L_X(D)`. In addition to the curve :math:`X`
and the divisor :math:`D`, you must also specify the evaluation
divisor :math:`E`.

Note that this section has not been updated since the wrapper
``riemann_roch_basis`` has been fixed. See above for how to
properly define a divisor for Singular's ``BrillNoether``
command.

Here's an example, one which computes a generator matrix of an
associated AG code. This time we use Singular's ``AGCode_L``
command.

::

    sage: singular.LIB('brnoeth.lib')
    sage: singular.eval("ring s = 2,(x,y),lp;")
    ''
    sage: print(singular.eval("list HC = Adj_div(x3+y2+y);"))
    Computing affine singular points ...
    Computing all points at infinity ...
    Computing affine singular places ...
    Computing singular places at infinity ...
    Computing non-singular places at infinity ...
    Adjunction divisor computed successfully
    <BLANKLINE>
    The genus of the curve is 1
    sage: print(singular.eval("list HC1 = NSplaces(1..2,HC);"))
    Computing non-singular affine places of degree 1 ...
    Computing non-singular affine places of degree 2 ...
    sage: print(singular.eval("HC = extcurve(2,HC1);"))
    Total number of rational places : NrRatPl = 9

We set the following to ``junk`` to discard the output::

    sage: junk = singular.eval("intvec G = 5;")      # the rational divisor G = 5*HC[3][1]
    sage: junk = singular.eval("def R = HC[1][2];")
    sage: singular.eval("setring R;")
    ''

The vector :math:`G` represents the divisor
"5 times the point at infinity".

.. index:: Riemann-Roch space

Next, we compute the Riemann-Roch space.

.. link

::

    sage: print(singular.eval("BrillNoether(G,HC);"))
    Forms of degree 3 :
    10
    <BLANKLINE>
    Vector basis successfully computed
    <BLANKLINE>
    [1]:
       _[1]=x
       _[2]=z
    [2]:
       _[1]=y
       _[2]=z
    [3]:
       _[1]=1
       _[2]=1
    [4]:
       _[1]=y2+yz
       _[2]=xz
    [5]:
       _[1]=y3+y2z
       _[2]=x2z

That was the basis of the Riemann-Roch space, where each pair of
functions represents the quotient (first function divided by second
function). Each of these basis elements get evaluated at certain
points to construct the generator matrix of the code. We next
construct the points.

.. skip

::

    sage: singular.eval("def R = HC[1][5];")
    '// ** redefining R **'
    sage: singular.eval("setring R;")
    ''
    sage: print(singular.eval("POINTS;"))
    [1]:
       [1]:
          0
       [2]:
          1
       [3]:
          0
    [2]:
       [1]:
          0
       [2]:
          1
       [3]:
          1
    [3]:
       [1]:
          0
       [2]:
          0
       [3]:
          1
    [4]:
       [1]:
          (a+1)
       [2]:
          (a)
       [3]:
          1
    ...

plus :math:`5` more, for a total of :math:`9` rational points
on the curve. We define our "evaluation divisor" :math:`D` using
a subset of these points (all but the first):

.. skip

::

    sage: singular.eval("def ER = HC[1][4];")
    ''
    sage: singular.eval("setring ER;")
    ''
    sage: # D = sum of the rational places no. 2..9 over F_4
    sage: singular.eval("intvec D = 2..9;")
    ''
    sage: # let us construct the corresponding evaluation AG code :
    sage: print(singular.eval("matrix C = AGcode_L(G,D,HC);"))
    Forms of degree 3 :
    10
    <BLANKLINE>
    Vector basis successfully computed
    <BLANKLINE>
    sage: # here is a linear code of type [8,5,> = 3] over F_4
    sage: print(singular.eval("print(C);"))
    0,0,(a+1),(a),  1,  1,    (a),  (a+1),
    1,0,(a),  (a+1),(a),(a+1),(a),  (a+1),
    1,1,1,    1,    1,  1,    1,    1,
    0,0,(a),  (a+1),1,  1,    (a+1),(a),
    0,0,1,    1,    (a),(a+1),(a+1),(a)

This is, finally, our desired generator matrix, where ``a``
represents a generator of the field extension of degree :math:`2`
over the base field :math:`GF(2)`.

Can this be "wrapped"?
