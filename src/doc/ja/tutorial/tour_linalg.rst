.. _section-linalg:

線形代数
==============

Sageには線形代数で常用されるツールが揃っていて，例えば行列の特性多項式の計算、階段形式、跡(トレース)、各種の分解などの操作が可能である．

行列の生成と積演算の手順は，簡単かつ自然なものだ:


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


Sageでは，行列 :math:`A` の核空間は「左」核空間，すなわち :math:`wA=0` を満足するベクトル :math:`w` が張る空間をさす．

..
   Note that in Sage, the kernel of a matrix :math:`A` is the
   "left kernel", i.e. the space of vectors :math:`w` such that
   :math:`wA=0`.

行列方程式もメソッド ``solve_right`` を使って簡単に解くことができる．
``A.solve_right(Y)`` と実行すれば， :math:`AX=Y` を満たす行列(またはベクトル) :math:`X` が得られる．


.. link

::

    sage: Y = vector([0, -4, -1])
    sage: X = A.solve_right(Y)
    sage: X
    (-2, 1, 0)
    sage: A * X   # 解のチェック...
    (0, -4, -1)


解がない場合は，Sageはエラーを返してくる:

.. skip

::

    sage: A.solve_right(w)
    Traceback (most recent call last):
    ...
    ValueError: matrix equation has no solutions


同様にして， :math:`XA=Y` を満足する :math:`X` を求めるには ``A.solve_left(Y)`` とすればよい．


Sageは固有値と固有ベクトルの計算もしてくれる:


::

    sage: A = matrix([[0, 4], [-1, 0]])
    sage: A.eigenvalues ()
    [-2*I, 2*I]
    sage: B = matrix([[1, 3], [3, 1]])
    sage: B.eigenvectors_left()
    [(4, [(1, 1)], 1), (-2, [(1, -1)], 1)]


( ``eigenvectors_left`` の出力は，三つ組タプル(固有値,固有ベクトル,多重度)のリストになっている．)

``QQ`` または ``RR`` 上の固有値と固有ベクトルはMaximaを使って計算することもできる( 後半の :ref:`section-maxima` 節を参照)．


:ref:`section-rings` 節で述べたように，行列の性質の中には，その行列がどんな環の上で定義されているかに影響を受けるものがある．
以下では， ``matrix`` コマンドの最初の引数を使って，Sageに生成すべき行列が整数の行列( ``ZZ`` の場合) なのか，有理数の行列(``QQ``)なのか，あるいは実数の行列(``RR``)なのかを指定している．


::

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


浮動小数点型の実数または複素数上で定義された行列の固有値と固有ベクトルを計算するためには，対象とする行列を，それぞれ ``RDF`` (Real Double Field)または ``CDF`` (Complex Double Field)上で定義しておかなければならない．
もし環を指定しないまま行列に浮動小数点型の実数あるいは複素数が使われる場合，その行列はデフォルトでそれぞれ ``RR`` あるいは ``CC`` 体上で定義される．
この場合，以下の演算があらゆる状況で実行可能になるとは限らない．


::

    sage: ARDF = matrix(RDF, [[1.2, 2], [2, 3]])
    sage: ARDF.eigenvalues()   # abs tol 1e-10
    [-0.09317121994613098, 4.293171219946131]
    sage: ACDF = matrix(CDF, [[1.2, I], [2, 3]])
    sage: ACDF.eigenvectors_right()   # abs tol 1e-10
    [(0.881845698329 - 0.820914065343*I, [(0.750560818381, -0.616145932705 + 0.238794153033*I)], 1),
    (3.31815430167 + 0.820914065343*I, [(0.145594698293 + 0.37566908585*I, 0.915245825866)], 1)]


行列の空間
-------------

有理数型の要素からなる `3 \times 3` 行列の空間
:math:`\text{Mat}_{3\times 3}(\QQ)` を生成してみよう:

::

    sage: M = MatrixSpace(QQ,3)
    sage: M
    Full MatrixSpace of 3 by 3 dense matrices over Rational Field


(3行4列の行列空間を生成したければ ``MatrixSpace(QQ,3,4)`` とする．
列数を省略するとデフォルトで行数に合わせられるから， ``MatrixSpace(QQ,3)`` は ``MatrixSpace(QQ,3,3)`` と同じ意味になる．)
行列の空間は標準基底を備えており:

.. link

::

    sage: B = M.basis()
    sage: len(B)
    9
    sage: B[0,1]
    [0 1 0]
    [0 0 0]
    [0 0 0]

``M`` の元の一つとして行列を生成してみよう．


.. link

::

    sage: A = M(range(9)); A
    [0 1 2]
    [3 4 5]
    [6 7 8]


ついで，その既約階段形式と核を計算する．

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

次に，有限体上で定義された行列による計算を実行してみる．


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

上に現れた行ベクトル系(rows)によって張られる `\GF{2}` の部分空間を作成する．


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


Sageは `S` の基底として， `S` の生成元行列の既約階段形式の非ゼロ行を使用している．


疎行列の線形代数
---------------------

SageではPID(単項イデアル整域)上の疎行列に関する線形代数を扱うことができる．

::

    sage: M = MatrixSpace(QQ, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()


Sageで使われている多重モジュラーアルゴリズムは，正方行列ではうまく働く(非正方行列ではいまひとつである):

::

    sage: M = MatrixSpace(QQ, 50, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()
    sage: M = MatrixSpace(GF(2), 20, 40, sparse=True)
    sage: A = M.random_element()
    sage: E = A.echelon_form()

Pythonでは，大文字小文字が区別されることに注意:

::

    sage: M = MatrixSpace(QQ, 10,10, Sparse=True)
    Traceback (most recent call last):
    ...
    TypeError: ...__init__() got an unexpected keyword argument 'Sparse'
