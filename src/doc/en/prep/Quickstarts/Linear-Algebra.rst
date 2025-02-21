.. -*- coding: utf-8 -*-

.. linkall

.. _prep-quickstart-linear-algebra:

Sage Quickstart for Linear Algebra
==================================

This `Sage <http://www.sagemath.org>`_ quickstart tutorial was developed
for the MAA PREP Workshop "Sage: Using Open\-Source Mathematics Software
with Undergraduates" (funding provided by NSF DUE 0817071).  It is
licensed under the Creative Commons Attribution\-ShareAlike 3.0 license
(`CC BY\-SA <http://creativecommons.org/licenses/by-sa/3.0/>`_).

Linear algebra underpins a lot of Sage's algorithms, so it is fast,
robust and comprehensive.  We've already seen some basic linear algebra,
including matrices, determinants, and the ``.rref()`` method for
row-reduced echelon form in the :doc:`Programming Tutorial
<../Programming>`, so the content here continues from there to some
extent.

Matrices and Vectors
--------------------

We can make a matrix easily by passing a list of the rows.  Don't forget
to use tab\-completion to see routines that are possible.

::

    sage: A = matrix([[1,2,3],[4,5,6]]); A
    [1 2 3]
    [4 5 6]

But there are lots of other ways to make matrices.  Each of these shows
what is assumed with different input; can you figure out how Sage
interprets them before you read the documentation which the command
``matrix?`` provides?

It's a good idea to get in the habit of telling Sage what ring to make
the matrix over.  Otherwise, Sage guesses based on the elements, so you
may not have a matrix over a field!  Here, we tell Sage to make the ring
over the rationals.

::

    sage: B = matrix(QQ, 3, 2, [1,2,3,4,5,6]); B
    [1 2]
    [3 4]
    [5 6]

::

    sage: C = matrix(QQ, 3, [1,2,3,4,5,6]); C
    [1 2]
    [3 4]
    [5 6]

::

    sage: D = matrix(CC, 20, range(400)); D
    20 x 20 dense matrix over Complex Field with 53 bits of precision (use the '.str()' method to see the entries)

Don't forget that when viewing this in the notebook, you can click to
the left of the matrix in order to cycle between "wrapped",
"unwrapped" and "hidden" modes of output.

::

    sage: print(D.str())
    [0.000000000000000  1.00000000000000  2.00000000000000  3.00000000000000  4.00000000000000  5.00000000000000  6.00000000000000  7.00000000000000  8.00000000000000  9.00000000000000  10.0000000000000  11.0000000000000  12.0000000000000  13.0000000000000  14.0000000000000  15.0000000000000  16.0000000000000  17.0000000000000  18.0000000000000  19.0000000000000]
    [ 20.0000000000000  21.0000000000000  22.0000000000000  23.0000000000000  24.0000000000000  25.0000000000000  26.0000000000000  27.0000000000000  28.0000000000000  29.0000000000000  30.0000000000000  31.0000000000000  32.0000000000000  33.0000000000000  34.0000000000000  35.0000000000000  36.0000000000000  37.0000000000000  38.0000000000000  39.0000000000000]
    [ 40.0000000000000  41.0000000000000  42.0000000000000  43.0000000000000  44.0000000000000  45.0000000000000  46.0000000000000  47.0000000000000  48.0000000000000  49.0000000000000  50.0000000000000  51.0000000000000  52.0000000000000  53.0000000000000  54.0000000000000  55.0000000000000  56.0000000000000  57.0000000000000  58.0000000000000  59.0000000000000]
    [ 60.0000000000000  61.0000000000000  62.0000000000000  63.0000000000000  64.0000000000000  65.0000000000000  66.0000000000000  67.0000000000000  68.0000000000000  69.0000000000000  70.0000000000000  71.0000000000000  72.0000000000000  73.0000000000000  74.0000000000000  75.0000000000000  76.0000000000000  77.0000000000000  78.0000000000000  79.0000000000000]
    [ 80.0000000000000  81.0000000000000  82.0000000000000  83.0000000000000  84.0000000000000  85.0000000000000  86.0000000000000  87.0000000000000  88.0000000000000  89.0000000000000  90.0000000000000  91.0000000000000  92.0000000000000  93.0000000000000  94.0000000000000  95.0000000000000  96.0000000000000  97.0000000000000  98.0000000000000  99.0000000000000]
    [ 100.000000000000  101.000000000000  102.000000000000  103.000000000000  104.000000000000  105.000000000000  106.000000000000  107.000000000000  108.000000000000  109.000000000000  110.000000000000  111.000000000000  112.000000000000  113.000000000000  114.000000000000  115.000000000000  116.000000000000  117.000000000000  118.000000000000  119.000000000000]
    [ 120.000000000000  121.000000000000  122.000000000000  123.000000000000  124.000000000000  125.000000000000  126.000000000000  127.000000000000  128.000000000000  129.000000000000  130.000000000000  131.000000000000  132.000000000000  133.000000000000  134.000000000000  135.000000000000  136.000000000000  137.000000000000  138.000000000000  139.000000000000]
    [ 140.000000000000  141.000000000000  142.000000000000  143.000000000000  144.000000000000  145.000000000000  146.000000000000  147.000000000000  148.000000000000  149.000000000000  150.000000000000  151.000000000000  152.000000000000  153.000000000000  154.000000000000  155.000000000000  156.000000000000  157.000000000000  158.000000000000  159.000000000000]
    [ 160.000000000000  161.000000000000  162.000000000000  163.000000000000  164.000000000000  165.000000000000  166.000000000000  167.000000000000  168.000000000000  169.000000000000  170.000000000000  171.000000000000  172.000000000000  173.000000000000  174.000000000000  175.000000000000  176.000000000000  177.000000000000  178.000000000000  179.000000000000]
    [ 180.000000000000  181.000000000000  182.000000000000  183.000000000000  184.000000000000  185.000000000000  186.000000000000  187.000000000000  188.000000000000  189.000000000000  190.000000000000  191.000000000000  192.000000000000  193.000000000000  194.000000000000  195.000000000000  196.000000000000  197.000000000000  198.000000000000  199.000000000000]
    [ 200.000000000000  201.000000000000  202.000000000000  203.000000000000  204.000000000000  205.000000000000  206.000000000000  207.000000000000  208.000000000000  209.000000000000  210.000000000000  211.000000000000  212.000000000000  213.000000000000  214.000000000000  215.000000000000  216.000000000000  217.000000000000  218.000000000000  219.000000000000]
    [ 220.000000000000  221.000000000000  222.000000000000  223.000000000000  224.000000000000  225.000000000000  226.000000000000  227.000000000000  228.000000000000  229.000000000000  230.000000000000  231.000000000000  232.000000000000  233.000000000000  234.000000000000  235.000000000000  236.000000000000  237.000000000000  238.000000000000  239.000000000000]
    [ 240.000000000000  241.000000000000  242.000000000000  243.000000000000  244.000000000000  245.000000000000  246.000000000000  247.000000000000  248.000000000000  249.000000000000  250.000000000000  251.000000000000  252.000000000000  253.000000000000  254.000000000000  255.000000000000  256.000000000000  257.000000000000  258.000000000000  259.000000000000]
    [ 260.000000000000  261.000000000000  262.000000000000  263.000000000000  264.000000000000  265.000000000000  266.000000000000  267.000000000000  268.000000000000  269.000000000000  270.000000000000  271.000000000000  272.000000000000  273.000000000000  274.000000000000  275.000000000000  276.000000000000  277.000000000000  278.000000000000  279.000000000000]
    [ 280.000000000000  281.000000000000  282.000000000000  283.000000000000  284.000000000000  285.000000000000  286.000000000000  287.000000000000  288.000000000000  289.000000000000  290.000000000000  291.000000000000  292.000000000000  293.000000000000  294.000000000000  295.000000000000  296.000000000000  297.000000000000  298.000000000000  299.000000000000]
    [ 300.000000000000  301.000000000000  302.000000000000  303.000000000000  304.000000000000  305.000000000000  306.000000000000  307.000000000000  308.000000000000  309.000000000000  310.000000000000  311.000000000000  312.000000000000  313.000000000000  314.000000000000  315.000000000000  316.000000000000  317.000000000000  318.000000000000  319.000000000000]
    [ 320.000000000000  321.000000000000  322.000000000000  323.000000000000  324.000000000000  325.000000000000  326.000000000000  327.000000000000  328.000000000000  329.000000000000  330.000000000000  331.000000000000  332.000000000000  333.000000000000  334.000000000000  335.000000000000  336.000000000000  337.000000000000  338.000000000000  339.000000000000]
    [ 340.000000000000  341.000000000000  342.000000000000  343.000000000000  344.000000000000  345.000000000000  346.000000000000  347.000000000000  348.000000000000  349.000000000000  350.000000000000  351.000000000000  352.000000000000  353.000000000000  354.000000000000  355.000000000000  356.000000000000  357.000000000000  358.000000000000  359.000000000000]
    [ 360.000000000000  361.000000000000  362.000000000000  363.000000000000  364.000000000000  365.000000000000  366.000000000000  367.000000000000  368.000000000000  369.000000000000  370.000000000000  371.000000000000  372.000000000000  373.000000000000  374.000000000000  375.000000000000  376.000000000000  377.000000000000  378.000000000000  379.000000000000]
    [ 380.000000000000  381.000000000000  382.000000000000  383.000000000000  384.000000000000  385.000000000000  386.000000000000  387.000000000000  388.000000000000  389.000000000000  390.000000000000  391.000000000000  392.000000000000  393.000000000000  394.000000000000  395.000000000000  396.000000000000  397.000000000000  398.000000000000  399.000000000000]

::

    sage: E = diagonal_matrix( [0..40,step=4] ); E
    [ 0  0  0  0  0  0  0  0  0  0  0]
    [ 0  4  0  0  0  0  0  0  0  0  0]
    [ 0  0  8  0  0  0  0  0  0  0  0]
    [ 0  0  0 12  0  0  0  0  0  0  0]
    [ 0  0  0  0 16  0  0  0  0  0  0]
    [ 0  0  0  0  0 20  0  0  0  0  0]
    [ 0  0  0  0  0  0 24  0  0  0  0]
    [ 0  0  0  0  0  0  0 28  0  0  0]
    [ 0  0  0  0  0  0  0  0 32  0  0]
    [ 0  0  0  0  0  0  0  0  0 36  0]
    [ 0  0  0  0  0  0  0  0  0  0 40]

::

    sage: column_matrix(QQ,[[1,2,3],[4,5,6],[7,8,9]])
    [1 4 7]
    [2 5 8]
    [3 6 9]

You can also combine matrices in different ways.

::

    sage: F1=matrix(QQ,2,2,[0,1,1,0])
    sage: F2=matrix(QQ,2,2,[1,2,3,4])
    sage: F3=matrix(QQ,1,2,[3,1])
    sage: block_matrix(2,2,[F1,F2,0,F3])
    [0 1|1 2]
    [1 0|3 4]
    [---+---]
    [0 0|3 1]

::

    sage: F1.augment(F2)
    [0 1 1 2]
    [1 0 3 4]

::

    sage: F1.stack(F2)
    [0 1]
    [1 0]
    [1 2]
    [3 4]

::

    sage: block_diagonal_matrix([F1,F2])
    [0 1|0 0]
    [1 0|0 0]
    [---+---]
    [0 0|1 2]
    [0 0|3 4]

Vectors are rows or columns, whatever you please, and Sage interprets
them as appropriate in multiplication contexts.

::

    sage: row = vector( (3, -1, 4))
    sage: col = vector( QQ, [4, 5] )
    sage: row; col
    (3, -1, 4)
    (4, 5)

::

    sage: F = matrix(QQ, 3, 2, range(6)); F
    [0 1]
    [2 3]
    [4 5]

::

    sage: F*col
    (5, 23, 41)

::

    sage: row*F
    (14, 20)

Although our "vectors" (especially over rings other than fields) might
be considered as elements of an appropriate free module, they basically behave as vectors
for our purposes.

::

    sage: ring_vec = vector(SR, [2, 12, -4, 9])
    sage: field_vec = vector( QQ, (2, 3, 14) )
    sage: ring_vec; field_vec
    (2, 12, -4, 9)
    (2, 3, 14)

::

    sage: type( ring_vec )
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category.element_class'>
    sage: type( field_vec )
    <class 'sage.modules.vector_rational_dense.Vector_rational_dense'>

Left\-Handed or Right\-handed?
-------------------------------

Sage "prefers" rows to columns.  For example, the ``kernel`` method
for a matrix `A` computes the left kernel -- the vector space of all
vectors `v` for which `v \cdot A = 0` -- and prints out the vectors as
the rows of a matrix.

::

    sage: G = matrix(QQ, 2, 3, [[1,2,3],[2,4,6]])
    sage: G.kernel()
    Vector space of degree 2 and dimension 1 over Rational Field
    Basis matrix:
    [   1 -1/2]

::

    sage: G.left_kernel()
    Vector space of degree 2 and dimension 1 over Rational Field
    Basis matrix:
    [   1 -1/2]

The ``right_kernel`` method computes the space of vectors `w` so that
`A \cdot w = 0`, of course.

Vector Spaces
--------------

Since Sage knows the kernel is a vector space, you can compute things
that make sense for a vector space.

::

    sage: V=G.right_kernel()
    sage: V
    Vector space of degree 3 and dimension 2 over Rational Field
    Basis matrix:
    [   1    0 -1/3]
    [   0    1 -2/3]

::

    sage: V.dimension()
    2

Here we compute the coordinate vector of :math:`(1,4,-3)` relative to
:math:`V`::

    sage: V.coordinate_vector([1,4,-3])
    (1, 4)

Here we get the basis matrix (note that the basis vectors are the *rows*
of the matrix)::

    sage: V.basis_matrix()
    [   1    0 -1/3]
    [   0    1 -2/3]

Or we can get the basis vectors explicitly as a list of vectors::

    sage: V.basis()
    [(1, 0, -1/3), (0, 1, -2/3)]

.. note::
   Kernels are **vector spaces** and bases are "\ **echelonized**\ "
   (canonicalized).

   This is why the ``ring`` for the matrix is important.  Compare the
   kernels above with the kernel using a matrix which is only defined over
   the integers.

   ::

       sage: G = matrix(ZZ,2, 3, [[1,2,3],[2,4,6]])
       sage: G.kernel()
       Free module of degree 2 and rank 1 over Integer Ring
       Echelon basis matrix:
       [ 2 -1]

Computations
-------------

Here are some more computations with matrices and vectors.

As you might expect, random matrices are random.

::

    sage: H = random_matrix(QQ, 5, 5, num_bound = 10, den_bound = 4)
    sage: H.det() # random
    15416
    sage: H.eigenvalues() # random
    [-10.08361801792048?, -2.682220984496031?, 4.739405672111427?, -1.320116668180795? - 10.88676412262347?*I, -1.320116668180795? + 10.88676412262347?*I]

According to the :doc:`Numerical analysis quickstart <NumAnalysis>`,
the question marks indicate that the actual
number is inside the interval found by incrementing and
decrementing the last digit of the printed number.  So 9.1? is a number
between 9.0 and 9.2.  Sage knows exactly what number this is (since it's
a root of a polynomial), but uses interval notation to print an
approximation for ease of use.

The ``eigenvectors_right`` command prints out a list of ``(eigenvalue,
[list of eigenvectors], algebraic multiplicity)`` tuples for each
eigenvalue.

::

    sage: H.eigenvectors_right() # random
    [(-10.08361801792048?, [(1, -0.3820692683963385?, -0.4659857618614747?, -0.1264082922197715?, -0.3548156445133095?)], 1), (-2.682220984496031?, [(1, -1.855347152382563?, -0.4203899923232704?, 0.004411201577480876?, -0.5050698736445243?)], 1), (4.739405672111427?, [(1, 0.3284800982819703?, 2.059182569319718?, -1.428547399599918?, 0.5455069936349178?)], 1), (-1.320116668180795? - 10.88676412262347?*I, [(1, 0.710831790589076? + 0.2646474741698805?*I, 0.4504038344112447? + 3.145667601780920?*I, 2.763061217778457? + 0.9994136057023008?*I, 3.092272491890536? - 2.105461094305392?*I)], 1), (-1.320116668180795? + 10.88676412262347?*I, [(1, 0.710831790589076? - 0.2646474741698805?*I, 0.4504038344112447? - 3.145667601780920?*I, 2.763061217778457? - 0.9994136057023008?*I, 3.092272491890536? + 2.105461094305392?*I)], 1)]

It may be more convenient to use the ``eigenmatrix_right`` command, which
gives a diagonal matrix of eigenvalues and a column matrix of
eigenvectors.

::

    sage: D,P=H.eigenmatrix_right()
    sage: P*D-H*P
    [0 0 0 0 0]
    [0 0 0 0 0]
    [0 0 0 0 0]
    [0 0 0 0 0]
    [0 0 0 0 0]

Matrix Solving
---------------

We can easily solve linear equations using the backslash, like in Matlab.

::

    sage: A = random_matrix(QQ, 3, algorithm='unimodular')
    sage: v = vector([2,3,1])
    sage: A,v  # random
    (
    [ 0 -1  1]
    [-1 -1 -1]
    [ 0  2  2], (2, 3, 1)
    )
    sage: x=A\v; x  # random
    (-7/2, -3/4, 5/4)
    sage: A*x  # random
    (2, 3, 1)

For *lots* more (concise) information, see the Sage `Linear Algebra
Quick Reference
<http://wiki.sagemath.org/quickref?action=AttachFile&do=get&target=quickref-linalg.pdf>`_.

