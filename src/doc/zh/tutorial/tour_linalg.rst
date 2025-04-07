.. _section-linalg:

线性代数
==============

Sage 提供了线性代数中的标准构造，例如矩阵的特征多项式、阶梯形、迹、分解等。

创建矩阵和进行矩阵乘法非常简单自然：

::

    sage: A = Matrix([[1,2,3],[3,2,1],[1,1,1]])
    sage: w = vector([1,1,-4])
    sage: w*A
    (0, 0, 0)
    sage: A*w
    (-9, 1, -2)
    sage: kernel(A)
    Free module of degree 3 and rank 1 over Integer Ring
    Echelon basis matrix:
    [ 1  1 -4]

请注意，在 Sage 中，矩阵 :math:`A` 的核是“左核”，即满足 :math:`wA=0` 的向量空间 :math:`w`。

求解矩阵方程非常简单，使用 ``solve_right`` 方法即可。
运行 ``A.solve_right(Y)`` 将返回一个矩阵（或向量） :math:`X`，使得 :math:`AX=Y`：

.. link

::

    sage: Y = vector([0, -4, -1])
    sage: X = A.solve_right(Y)
    sage: X
    (-2, 1, 0)
    sage: A * X   # checking our answer...
    (0, -4, -1)

倘若无解，Sage 会返回错误：

.. skip

::

    sage: A.solve_right(w)
    Traceback (most recent call last):
    ...
    ValueError: matrix equation has no solutions

同理，可以使用 ``A.solve_left(Y)`` 来求解方程 ::math:`XA=Y` 中的 :math:`X`。

Sage 还可以计算特征值和特征向量::

    sage: A = matrix([[0, 4], [-1, 0]])
    sage: A.eigenvalues ()
    [-2*I, 2*I]
    sage: B = matrix([[1, 3], [3, 1]])
    sage: B.eigenvectors_left()
    [(4, [(1, 1)], 1), (-2, [(1, -1)], 1)]

（``eigenvectors_left`` 的输出格式是一个包含三元组（特征值、特征向量、重数）的列表。）
特征值和特征向量可以通过 Maxima 在有理数域 ``QQ`` 或实数域 ``RR`` 上计算（见下文的 :ref:`section-maxima`）。

如 :ref:`section-rings` 所述，矩阵定义的环会影响其某些性质。
在下面的示例中，``matrix`` 命令的第一个参数告诉 Sage 将矩阵视为整数矩阵（``ZZ``）、有理数矩阵（``QQ``）或实数矩阵（``RR``） ::

    sage: AZ = matrix(ZZ, [[2,0], [0,1]])
    sage: AQ = matrix(QQ, [[2,0], [0,1]])
    sage: AR = matrix(RR, [[2,0], [0,1]])
    sage: AZ.echelon_form()
    [2 0]
    [0 1]
    sage: AQ.echelon_form()
    [1 0]
    [0 1]
    sage: AR.echelon_form()
    [ 1.00000000000000 0.000000000000000]
    [0.000000000000000  1.00000000000000]

如果要计算浮点实数或复数矩阵的特征值和特征向量，矩阵应分别定义在 ``RDF`` （实双精度域）或 ``CDF`` （复双精度域）上。
如果没有指定环并且使用浮点实数或复数，则默认情况下矩阵定义在 ``RR`` 或 ``CC`` 域上，这些域不支持所有情况的这些计算::

    sage: ARDF = matrix(RDF, [[1.2, 2], [2, 3]])
    sage: ARDF.eigenvalues()  # rel tol 8e-16
    [-0.09317121994613098, 4.293171219946131]
    sage: ACDF = matrix(CDF, [[1.2, I], [2, 3]])
    sage: ACDF.eigenvectors_right()  # rel tol 3e-15
    [(0.8818456983293743 - 0.8209140653434135*I, [(0.7505608183809549, -0.616145932704589 + 0.2387941530333261*I)], 1),
    (3.3181543016706256 + 0.8209140653434133*I, [(0.14559469829270957 + 0.3756690858502104*I, 0.9152458258662108)], 1)]

矩阵空间
-------------

我们创建了一个定义在有理数域 :math:`\QQ` 上的 `3 \times 3` 矩阵空间 :math:`\text{Mat}_{3\times 3}(\QQ)`::

    sage: M = MatrixSpace(QQ,3)
    sage: M
    Full MatrixSpace of 3 by 3 dense matrices over Rational Field

（要创建一个 `3 \times 4` 矩阵空间，可以使用 ``MatrixSpace(QQ,3,4)``。
如果省略列数，则默认为行数，因此 ``MatrixSpace(QQ,3)`` 与 ``MatrixSpace(QQ,3,3)`` 意义相同。）
矩阵空间有其规范基：

.. link

::

    sage: B = M.basis()
    sage: len(B)
    9
    sage: B[0,1]
    [0 1 0]
    [0 0 0]
    [0 0 0]

我们创建一个矩阵作为 ``M`` 的元素。

.. link

::

    sage: A = M(range(9)); A
    [0 1 2]
    [3 4 5]
    [6 7 8]

接下来我们计算其简化行阶梯形和核。

.. link

::

    sage: A.echelon_form()
    [ 1  0 -1]
    [ 0  1  2]
    [ 0  0  0]
    sage: A.kernel()
    Vector space of degree 3 and dimension 1 over Rational Field
    Basis matrix:
    [ 1 -2  1]

接着我们来演示在有限域上定义的矩阵的计算：

::

    sage: M = MatrixSpace(GF(2),4,8)
    sage: A = M([1,1,0,0, 1,1,1,1, 0,1,0,0, 1,0,1,1,
    ....:        0,0,1,0, 1,1,0,1, 0,0,1,1, 1,1,1,0])
    sage: A
    [1 1 0 0 1 1 1 1]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 1 1 1 1 1 0]
    sage: rows = A.rows()
    sage: A.columns()
    [(1, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 1), (0, 0, 0, 1),
     (1, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0)]
    sage: rows
    [(1, 1, 0, 0, 1, 1, 1, 1), (0, 1, 0, 0, 1, 0, 1, 1),
     (0, 0, 1, 0, 1, 1, 0, 1), (0, 0, 1, 1, 1, 1, 1, 0)]

我们创建一个在有限域 `\GF{2}` 上由上述行生成的子空间。

.. link

::

    sage: V = VectorSpace(GF(2),8)
    sage: S = V.subspace(rows)
    sage: S
    Vector space of degree 8 and dimension 4 over Finite Field of size 2
    Basis matrix:
    [1 0 0 0 0 1 0 0]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 0 1 0 0 1 1]
    sage: A.echelon_form()
    [1 0 0 0 0 1 0 0]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 0 1 0 0 1 1]

Sage 使用的 `S` 的基是通过生成矩阵的简化行阶梯形的非零行获得的。

稀疏线性代数
---------------------

Sage 支持在主理想域 (PIDs) 上的稀疏线性代数。

::

    sage: M = MatrixSpace(QQ, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()

Sage 中的多模算法适用于方阵（但不适用于非方阵）：

::

    sage: M = MatrixSpace(QQ, 50, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()
    sage: M = MatrixSpace(GF(2), 20, 40, sparse=True)
    sage: A = M.random_element()
    sage: E = A.echelon_form()

请注意，Python 是区分大小写的：

::

    sage: M = MatrixSpace(QQ, 10,10, Sparse=True)
    Traceback (most recent call last):
    ...
    TypeError: ...__init__() got an unexpected keyword argument 'Sparse'...
