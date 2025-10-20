.. _chapter-linear_algebra:

********
线性代数
********

.. index:
   pair: vector space; basis
   pair: vector space; subspace

.. _section-vector_space:

向量空间
========

``VectorSpace`` 命令用于创建一个向量空间类，从中可以创建一个子空间。
注意，Sage 计算的基是“行最简”的。

::

    sage: V = VectorSpace(GF(2),8)
    sage: S = V.subspace([V([1,1,0,0,0,0,0,0]),V([1,0,0,0,0,1,1,0])])
    sage: S.basis()
    [(1, 0, 0, 0, 0, 1, 1, 0), (0, 1, 0, 0, 0, 1, 1, 0)]
    sage: S.dimension()
    2


.. index:
   pair: matrix; powers

.. _section-matrixpower:

矩阵幂
======

如何在 Sage 中计算矩阵幂？语法如下例所示。

::

    sage: R = IntegerModRing(51)
    sage: M = MatrixSpace(R,3,3)
    sage: A = M([1,2,3, 4,5,6, 7,8,9])
    sage: A^1000*A^1007
    <BLANKLINE>
    [ 3  3  3]
    [18  0 33]
    [33 48 12]
    sage: A^2007
    <BLANKLINE>
    [ 3  3  3]
    [18  0 33]
    [33 48 12]

.. index:
   pair: matrix; kernel
   single: kernel; nullspace

.. _section-kernel:

核
==

通过对矩阵对象应用 ``kernel`` 方法来计算核。语法如下例所示。

::

    sage: M = MatrixSpace(IntegerRing(),4,2)(range(8))
    sage: M.kernel()
    Free module of degree 4 and rank 2 over Integer Ring
    Echelon basis matrix:
    [ 1  0 -3  2]
    [ 0  1 -2  1]

`\QQ` 上的一维核：

::

    sage: A = MatrixSpace(RationalField(),3)(range(9))
    sage: A.kernel()
    Vector space of degree 3 and dimension 1 over Rational Field
    Basis matrix:
    [ 1 -2  1]

平凡核：

::

    sage: A = MatrixSpace(RationalField(),2)([1,2,3,4])
    sage: A.kernel()
    Vector space of degree 2 and dimension 0 over Rational Field
    Basis matrix:
    []
    sage: M = MatrixSpace(RationalField(),0,2)(0)
    sage: M
    []
    sage: M.kernel()
    Vector space of degree 0 and dimension 0 over Rational Field
    Basis matrix:
    []
    sage: M = MatrixSpace(RationalField(),2,0)(0)
    sage: M.kernel()
    Vector space of degree 2 and dimension 2 over Rational Field
    Basis matrix:
    [1 0]
    [0 1]

零矩阵的核：

::

    sage: A = MatrixSpace(RationalField(),2)(0)
    sage: A.kernel()
    Vector space of degree 2 and dimension 2 over Rational Field
    Basis matrix:
    [1 0]
    [0 1]

非方阵的核：

::

    sage: A = MatrixSpace(RationalField(),3,2)(range(6))
    sage: A.kernel()
    Vector space of degree 3 and dimension 1 over Rational Field
    Basis matrix:
    [ 1 -2  1]

分圆域上矩阵的二维核：

::

    sage: K = CyclotomicField(12); a = K.gen()
    sage: M = MatrixSpace(K,4,2)([1,-1, 0,-2, 0,-a^2-1, 0,a^2-1])
    sage: M
    [             1            -1]
    [             0            -2]
    [             0 -zeta12^2 - 1]
    [             0  zeta12^2 - 1]
    sage: M.kernel()
    Vector space of degree 4 and dimension 2 over Cyclotomic Field of order 12
     and degree 4
    Basis matrix:
    [               0                1                0     -2*zeta12^2]
    [               0                0                1 -2*zeta12^2 + 1]

复杂基域上的非平凡核。

::

    sage: K = FractionField(PolynomialRing(RationalField(),2,'x'))
    sage: M = MatrixSpace(K, 2)([[K.gen(1),K.gen(0)], [K.gen(1), K.gen(0)]])
    sage: M
    [x1 x0]
    [x1 x0]
    sage: M.kernel()
    Vector space of degree 2 and dimension 1 over Fraction Field of Multivariate
    Polynomial Ring in x0, x1 over Rational Field
    Basis matrix:
     [ 1 -1]

.. index:: Smith normal form, Hermite normal form, Frobenius normal form, rational canonical form

其他一些适用于整数矩阵的方法包括 ``elementary_divisors``, ``smith_form``
（用于 Smith 标准型）, ``echelon_form`` 用于 Hermite 标准型，
``frobenius`` 用于 Frobenius 标准型（有理规范型）。


有许多域（例如：`\QQ`）或有限域上的矩阵方法：
``row_span``, ``nullity``,
``transpose``, ``swap_rows``, ``matrix_from_columns``,
``matrix_from_rows`` 等等。

请参阅文件 ``matrix.py`` 了解更多详情。

.. index:: eigenvalues, eigenvectors

.. _section-eigen:

特征向量和特征值
================

如何使用 Sage 计算特征值和特征向量？

Sage 提供了一系列完整的函数来计算特征值和左右特征向量以及特征子空间。
如果我们的矩阵为 `A`，那么 ``eigenmatrix_right``
（相对的为 ``eightmatrix_left``）命令会给出矩阵 `D` 和 `P`，
使得 `AP=PD`（相对的为 `PA=DP`）。

::

    sage: A = matrix(QQ, [[1,1,0],[0,2,0],[0,0,3]])
    sage: A
    [1 1 0]
    [0 2 0]
    [0 0 3]
    sage: A.eigenvalues()
    [3, 2, 1]
    sage: A.eigenvectors_right()
    [(3, [(0, 0, 1)], 1), (2, [(1, 1, 0)], 1), (1, [(1, 0, 0)], 1)]

    sage: A.eigenspaces_right()
    [(3,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [0 0 1]),
     (2,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [1 1 0]),
     (1,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [1 0 0])]

    sage: D, P = A.eigenmatrix_right()
    sage: D
    [3 0 0]
    [0 2 0]
    [0 0 1]
    sage: P
    [0 1 1]
    [0 1 0]
    [1 0 0]
    sage: A*P == P*D
    True

对于矩阵基环分数域之外的特征值，可以选择在实现域的代数闭包（例如代数数 ``QQbar``）时，
输出所有特征子空间。也可以为特征多项式的每个不可约因子请求独立的特征子空间，
因为其他特征子空间可以通过伽罗瓦共轭形成。下列矩阵的特征值为 $\pm\sqrt{3}$，
我们展示了每种可能的输出。

此外，目前 Sage 尚未实现多重精度数值特征值和特征向量，因此在 ``CC`` 或 ``RR`` 上
调用特征函数可能会给出不准确且无意义的结果（还会打印警告）。具有浮点项的矩阵
（在 ``CDF`` 和 ``RDF`` 上）的特征值和特征向量可以通过 "eigenmatrix" 命令获得。  ::

    sage: MS = MatrixSpace(QQ, 2, 2)
    sage: A = MS([1,-4,1, -1])
    sage: A.eigenspaces_left(format='all')
    [(-1.732050807568878?*I,
      Vector space of degree 2 and dimension 1 over Algebraic Field
      User basis matrix:
      [                        1 -1 - 1.732050807568878?*I]),
     (1.732050807568878?*I,
      Vector space of degree 2 and dimension 1 over Algebraic Field
      User basis matrix:
      [                        1 -1 + 1.732050807568878?*I])]
    sage: A.eigenspaces_left(format='galois')
    [(a0,
      Vector space of degree 2 and dimension 1 over Number Field in a0 with defining polynomial x^2 + 3
      User basis matrix:
      [     1 a0 - 1])]

另一种方法是通过接口调用 Maxima：

::

    sage: A = maxima("matrix ([1, -4], [1, -1])")
    sage: eig = A.eigenvectors()
    sage: eig.sage()
    [[[-I*sqrt(3), I*sqrt(3)], [1, 1]], [[[1, 1/4*I*sqrt(3) + 1/4]], [[1, -1/4*I*sqrt(3) + 1/4]]]]

这告诉我们 `\vec{v}_1 = [1,(\sqrt{3}i + 1)/4]` 是
`\lambda_1 = - \sqrt{3}i` （重数为 1）的特征向量，
而 `\vec{v}_2 = [1,(-\sqrt{3}i + 1)/4]` 是
`\lambda_2 =  \sqrt{3}i` （重数也为 1）的特征向量。

以下是另外两个例子：

::

    sage: A = maxima("matrix ([11, 0, 0], [1, 11, 0], [1, 3, 2])")
    sage: A.eigenvectors()
    [[[2,11],[1,2]],[[[0,0,1]],[[0,1,1/3]]]]
    sage: A = maxima("matrix ([-1, 0, 0], [1, -1, 0], [1, 3, 2])")
    sage: A.eigenvectors()
    [[[-1,2],[2,1]],[[[0,1,-1]],[[0,0,1]]]]

警告：请注意输出的顺序是相反的，尽管矩阵几乎相同。

最后，你还可以使用 Sage 的 GAP 接口来计算“有理”特征值和特征向量：

::

    sage: A = libgap([[1,2,3],[4,5,6],[7,8,9]]); A
    [ [ 1, 2, 3 ], [ 4, 5, 6 ], [ 7, 8, 9 ] ]
    sage: libgap(QQ).Eigenvectors(A)
    [ [ 1, -2, 1 ] ]
    sage: libgap(QQ).Eigenvalues(A)
    [ 0 ]

.. _section-rref:

行化简
======

矩阵的行最简阶梯形式按以下示例计算。

::

    sage: M = MatrixSpace(RationalField(),2,3)
    sage: A = M([1,2,3, 4,5,6])
    sage: A
    [1 2 3]
    [4 5 6]
    sage: A.parent()
    Full MatrixSpace of 2 by 3 dense matrices over Rational Field
    sage: A[0,2] = 389
    sage: A
    [  1   2 389]
    [  4   5   6]
    sage: A.echelon_form()
    [      1       0 -1933/3]
    [      0       1  1550/3]

.. index::
   pair: matrix; characteristic polynomial

.. _section-characteristic:

特征多项式
==========

特征多项式是一个适用于方阵的 Sage 方法。

首先是 `\ZZ` 上的矩阵：

::

    sage: A = MatrixSpace(IntegerRing(),2)( [[1,2], [3,4]] )
    sage: f = A.charpoly()
    sage: f
    x^2 - 5*x - 2
    sage: f.parent()
    Univariate Polynomial Ring in x over Integer Ring

我们计算一个定义在多项式环 `\ZZ[a]` 上的矩阵的特征多项式：

::

    sage: R = PolynomialRing(IntegerRing(),'a'); a = R.gen()
    sage: M = MatrixSpace(R,2)([[a,1], [a,a+1]])
    sage: M
    [    a     1]
    [    a a + 1]
    sage: f = M.charpoly()
    sage: f
    x^2 + (-2*a - 1)*x + a^2
    sage: f.parent()
    Univariate Polynomial Ring in x over Univariate Polynomial Ring in a over
    Integer Ring

    sage: M.trace()
    2*a + 1
    sage: M.determinant()
    a^2

我们计算一个定义在多元多项式环 `\ZZ[u,v]` 上的矩阵的特征多项式：

::

    sage: R.<u,v> = PolynomialRing(ZZ,2)
    sage: A = MatrixSpace(R,2)([u,v,u^2,v^2])
    sage: f = A.charpoly(); f
    x^2 + (-v^2 - u)*x - u^2*v + u*v^2

区分变量有点困难。为了解决这个问题，我们可能需要重命名不定变量 "Z"，我们可以轻松地执行如下操作：

.. link

::

    sage: f = A.charpoly('Z'); f
    Z^2 + (-v^2 - u)*Z - u^2*v + u*v^2

.. index::
   pair: solve; linear equations

求解线性方程组
==============

使用 Maxima，可以轻松求解线性方程：

::

    sage: var('a,b,c')
    (a, b, c)
    sage: eqn = [a+b*c==1, b-a*c==0, a+b==5]
    sage: s = solve(eqn, a,b,c); s
    [[a == -1/4*I*sqrt(79) + 11/4, b == 1/4*I*sqrt(79) + 9/4, c == 1/10*I*sqrt(79) + 1/10], [a == 1/4*I*sqrt(79) + 11/4, b == -1/4*I*sqrt(79) + 9/4, c == -1/10*I*sqrt(79) + 1/10]]

你甚至可以用 LaTeX 很好地排版方程的解：

::

    sage.: print(latex(s))
    ...

要通过 xdvi 在屏幕上显示上述内容，请输入 ``view(s)``。

你还可以通过 ``solve`` 命令符号化求解线性方程::

    sage: var('x,y,z,a')
    (x, y, z, a)
    sage: eqns = [x + z == y, 2*a*x - y == 2*a^2, y - 2*z == 2]
    sage: solve(eqns, x, y, z)
    [[x == a + 1, y == 2*a, z == a - 1]]

这是一个数值化 Numpy 示例::

    sage: from numpy import arange, eye, linalg
    sage: A = eye(10)       ##   the 10x10 identity matrix
    sage: b = arange(1,11)
    sage: x = linalg.solve(A,b)

另一种数值化求解方程组的方法是使用 Sage 的 Octave 接口::

    sage: M33 = MatrixSpace(QQ,3,3)
    sage: A   = M33([1,2,3,4,5,6,7,8,0])
    sage: V3  = VectorSpace(QQ,3)
    sage: b   = V3([1,2,3])
    sage: octave.solve_linear_system(A,b)    # optional - octave
    [-0.333333, 0.666667, 0]
