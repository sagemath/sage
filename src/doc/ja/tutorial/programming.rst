=================
 プログラミング
=================

.. _section-loadattach:

Sageファイルの読み込みと結合
==============================

ここでは，独立したファイルに保存したプログラムをSageに読み込む方法を解説する．
まず，ファイル ``example.sage`` に以下のプログラムを保存しておこう:

.. skip

::

    print("Hello World")
    print(2^3)

``example.sage`` ファイルを読み込んで実行するには， ``load`` コマンド を使う．

.. skip

::

    sage: load("example.sage")
    Hello World
    8

実行中のセッションにSageファイルを結合するには， ``attach`` コマンドが使える:

.. skip

::

    sage: attach("example.sage")
    Hello World
    8

こうして ``attach`` したファイル ``example.sage`` に変更を加えてからSageで空行を入力すると( ``return`` を押す)，Sageは自動的に ``example.sage`` の内容を読み込んで実行する．

``attach`` は結合したファイルに変更が生じると自動的に読込んで実行してくれるので，デバッグの時にはとりわけ便利なコマンドになる．
これに対し ``load`` の方は，コマンド実行時に一度だけファイルを読込んで実行するだけだ．

``example.sage`` を読込むとSageは内容をPythonプログラムへと変換し，これをPythonインタプリタが実行する．
このPythonへの変換は最小限に留められていて，整数リテラルを ``Integer()`` で，浮動小数点数リテラルを ``RealNumber()`` でラップし， ``^`` を ``**`` で，また ``R.2`` を ``R.gen(2)`` で置換するなどといった程度である．
変換された ``example.sage`` のコードは， ``example.sage`` と同じディレクトリに ``example.sage.py`` という名前のファイルとして保存されている．
この ``example.sage.py`` に入っているコードは:

::

    print("Hello World")
    print(Integer(2)**Integer(3))

たしかに整数リテラルはラップされ， ``^`` は ``**`` に置換されている．(Pythonでは ``^`` は「排他的論理和」， ``**`` は「べき乗」を意味する．)

こうした前処理は ``sage/misc/interpreter.py`` として実装されている．

コードブロックが改行で開始されているなら，インデントされた複数行のコードをSageにペーストすることができる(ファイルの読込みについてはブロック開始時の改行も不要だ)．
しかしSage上にコードを取り込むには，そのコードをファイルに保存してから，先に説明したように ``attach`` するのが一番安全だ．


.. _section-compile:

実行形式の作成
===============

数学的演算処理では実行速度がとりわけ重要になる．
Pythonは使い勝手のよい非常に高水準の言語ではあるが，計算の種類によってはスタティックに型付けされたコンパイラ言語で実行した方が処理速度が数桁も速くなる場合がある．
Sageにも，もし完全にPythonのみで実装していたら遅くて使いものにならなかったはずの機能がある．
そうした状況を切り抜けるために，SageはCythonと呼ばれるコンパイラ言語版Pythonを利用している([Cyt]_ と [Pyr]_)．
Cythonは，PythonとCの両方に似ていて，リスト内包表記、条件制御、それに ``+=`` などを含むPythonの構文の大半を使うことができるし，Pythonで書いておいたモジュールをインポートすることも可能だ．
加えて，C言語の変数を自由に宣言し，かつ好きなCライブラリをコード内から直接呼び出すことができる．
Cythonで書いたコードはCに変換され，Cコンパイラでコンパイルされることになる．

Sageコードの実行形式を作成するには，ソースファイルの拡張子を(``.sage`` ではなく) ``.spyx`` とする．
コマンドラインインターフェイスを使っている場合，実行形式をインタプリタ用コードのときと全く同じやり方で ``attach`` あるいは ``load`` することができる(今のところ，ノートブックインターフェイスではCythonで書いたコードの読込みと結合はできない)．
実際のコンパイル処理は，ユーザーが特に何も指定しなくとも裏で実行されている．
作成された実行形式の共有ライブラリは ``$HOME/.sage/temp/hostname/pid/spyx`` に格納されるが，Sageの終了時に消去される．


Sageはspyxファイルに対しては前処理をしない．
だから例えばspyxファイル中の ``1/3`` は有理数 :math:`1/3` ではなく0と評価される．Sageライブラリに含まれる関数 ``foo`` をspyxファイルで使用したければ， ``sage.all`` をインポートしてから ``sage.all.foo`` として呼び出す必要がある．

::

    import sage.all
    def foo(n):
        return sage.all.factorial(n)


他ファイル中のC関数を使う
-------------------------

別途 \*.cファイルで定義されたC関数を使うのも簡単だ．
実例を見てみよう．
まず，同じディレクトリにファイル ``test.c`` と ``test.spyx`` を作り，内容はそれぞれ以下の通りであるとする:


純粋にCで書かれたプログラムは ``test.c`` :

::

    int add_one(int n) {
      return n + 1;
    }

Cythonプログラムが入っているのが ``test.spyx``:

::

    cdef extern from "test.c":
        int add_one(int n)

    def test(n):
        return add_one(n)

すると，次の手順でCプログラムを取り込むことができる:

.. skip

::

    sage: attach("test.spyx")
    Compiling (...)/test.spyx...
    sage: test(10)
    11

Cythonソースファイルから生成されたC言語コードをコンパイルするために，さらにライブラリ ``foo`` が必要な場合は，元のCythonソースに ``clib foo`` という行を加える．
同様に、他にもCソースファイル ``bar`` が必要ならばCythonソースに取り込みを宣言する行 ``cfile bar`` を加えてコンパイルすればよい．


.. _section-standalone:

スタンドアロンPython/Sageスクリプト
====================================

以下のスタンドアロン型Sageスクリプトは，整数や多項式などを因数分解する:

::

    #!/usr/bin/env sage

    import sys

    if len(sys.argv) != 2:
        print("Usage: %s <n>" % sys.argv[0])
        print("Outputs the prime factorization of n.")
        sys.exit(1)

    print(factor(sage_eval(sys.argv[1])))

このスクリプトを実行するには， ``SAGE_ROOT`` をPATHに含めておかなければならない．
スクリプト名を ``factor`` とすると，実行は以下のような具合になる:


::

    bash $ ./factor 2006
    2 * 17 * 59


データ型
=========

Sageに現れるオブジェクトには，全て明確に定義されたデータ型が割り当てられている．
Pythonは豊富な組み込み型を備えているが，それをさらに多彩に拡張しているのがSageのライブラリだ．
Pythonの組み込み型としては，string(文字列)，list(リスト)，タプル(tuple)，int(整数)，float(浮動小数点数)などがある．
実際に型を表示してみると:


::

    sage: s = "sage"; type(s)
    <... 'str'>
    sage: s = 'sage'; type(s)      # シングルあるいはダブル クォーテーションのどちらも使える
    <... 'str'>
    sage: s = [1,2,3,4]; type(s)
    <... 'list'>
    sage: s = (1,2,3,4); type(s)
    <... 'tuple'>
    sage: s = int(2006); type(s)
    <... 'int'>
    sage: s = float(2006); type(s)
    <... 'float'>

Sageでは，さらに多様な型が加わる．
その一例がベクトル空間である:

::

    sage: V = VectorSpace(QQ, 1000000); V
    Vector space of dimension 1000000 over Rational Field
    sage: type(V)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category'>

この ``V`` に適用できるのは，あらかじめ定められた特定の関数に限られる．
他の数学ソフトウェアでは，そうした関数の呼び出しに「関数型記法」 ``foo(V,...)`` が用いられているようだ．
これに対し，Sageでは ``V`` の型(クラス)に付属する関数群が定められていて，JAVAやC++に見られるようなオブジェクト指向型の構文 ``V.foo(...)`` で関数呼び出しが行なわれる．
オブジェクト指向の世界では，関数が莫大な数になってもグローバルな名前空間を混乱なく運用することが可能になる．
たとえ機能の異なる関数群が同じ ``foo`` という名前を持っていたとしても，型チェックや場合分け抜きで引数の型に応じた適切な関数が自動的に呼び出されるのである．
さらに，ある関数の名前を他の意味で使い回しても関数そのものは使い続けることができる(例えば ``zeta`` を何かの変数名として使った後でも，リーマンのゼータ関数の0.5における値を求めるには ``s=.5; s.zeta()`` と入力すれば足りる)．


::

    sage: zeta = -1
    sage: s=.5; s.zeta()
    -1.46035450880959


ごく慣用化している場合については，通常の関数型記法を使うこともできる．
これは便利であると同時に，数学表現にはそもそもオブジェクト指向型記法になじまないものもあるためだ．
ここで少し例を見てみよう．


::

    sage: n = 2; n.sqrt()
    sqrt(2)
    sage: sqrt(2)
    sqrt(2)
    sage: V = VectorSpace(QQ,2)
    sage: V.basis()
        [
        (1, 0),
        (0, 1)
        ]
    sage: basis(V)
        [
        (1, 0),
        (0, 1)
        ]
    sage: M = MatrixSpace(GF(7), 2); M
    Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
    sage: A = M([1,2,3,4]); A
    [1 2]
    [3 4]
    sage: A.charpoly('x')
    x^2 + 2*x + 5
    sage: charpoly(A, 'x')
    x^2 + 2*x + 5

:math:`A` のメンバ関数を全て表示するには，タブ補完入力を利用すればよい．
:ref:`section-tabcompletion` 節で説明したように，これは ``A.`` と入力してから，キーボードの ``[tab]`` キーを押すだけのことだ．


リスト，タプル，シーケンス
===========================

リスト型には，任意の型の要素を格納することができる．
(大半のコンピュータ代数システムとは違って)CやC++などと同じように，Sageでもリストの要素番号は :math:`0` から始まる:


::

    sage: v = [2, 3, 5, 'x', SymmetricGroup(3)]; v
    [2, 3, 5, 'x', Symmetric group of order 3! as a permutation group]
    sage: type(v)
    <... 'list'>
    sage: v[0]
    2
    sage: v[2]
    5

リストの要素番号は，Pythonのint型でなくとも平気だ．
SageのIntegerクラスが使えるのは言うまでもない(Rationalクラスを含めて、 ``__index__``  メソッドが有効なクラスであれば何でも使える)．

.. ..
..    (When indexing into a list, it is OK if the index is
..    not a Python int!)
..    A Sage Integer (or Rational, or anything with an ``__index__`` method)
..    will work just fine.

::

    sage: v = [1,2,3]
    sage: v[2]
    3
    sage: n = 2      # SAGEの整数
    sage: v[n]       # 全く問題なし
    3
    sage: v[int(n)]  # これも大丈夫
    3

``range`` 関数は，Pythonのint型からなるリストを生成する(SageのIntegerではないことに注意):

::

    sage: list(range(1, 15))
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

この ``range`` が便利なのは，リスト内包表記を使ってリストを生成する場合だ:


::

    sage: L = [factor(n) for n in range(1, 15)]
    sage: L
    [1, 2, 3, 2^2, 5, 2 * 3, 7, 2^3, 3^2, 2 * 5, 11, 2^2 * 3, 13, 2 * 7]
    sage: L[12]
    13
    sage: type(L[12])
    <class 'sage.structure.factorization_integer.IntegerFactorization'>
    sage: [factor(n) for n in range(1, 15) if is_odd(n)]
    [1, 3, 5, 7, 3^2, 11, 13]

以上のようなリスト内包表記を使ったリスト生成については， [PyT]_ に詳しい．

.. ..
..    For more about how to create lists using list comprehensions, see
..    [PyT]_.

とりわけ使い勝手が良いのが，リストのスライシングだ．
リスト ``L`` のスライシング ``L[m:n]`` は， :math:`m` 番目の要素に始まり :math:`n-1` 番目の要素で終わる部分リストを返す．
以下に例を示そう:

::

    sage: L = [factor(n) for n in range(1, 20)]
    sage: L[4:9]
    [5, 2 * 3, 7, 2^3, 3^2]
    sage: L[:4]
    [1, 2, 3, 2^2]
    sage: L[14:4]
    []
    sage: L[14:]
    [3 * 5, 2^4, 17, 2 * 3^2, 19]

タプルはリストに似ているが，これがいったん生成された後は変更できない不変性(immutable)オブジェクトである点で異なる．


::

    sage: v = (1,2,3,4); v
    (1, 2, 3, 4)
    sage: type(v)
    <... 'tuple'>
    sage: v[1] = 5
    Traceback (most recent call last):
    ...
    TypeError: 'tuple' object does not support item assignment

Sageで使われる第三のリスト類似データ型が，シーケンスである．
リストやタプルと違って，シーケンスはPython本体の組み込み型ではない．
デフォルトではシーケンス型は可変だが，以下の例で見るように ``Sequence`` クラスのメソッド ``set_immutable`` を使って不変性を与えることができる．
あるシーケンスの全要素は，シーケンス・ユニバースと呼ばれる共通のペアレントを持つ．


::

    sage: v = Sequence([1,2,3,4/5])
    sage: v
    [1, 2, 3, 4/5]
    sage: type(v)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: type(v[1])
    <class 'sage.rings.rational.Rational'>
    sage: v.universe()
    Rational Field
    sage: v.is_immutable()
    False
    sage: v.set_immutable()
    sage: v[0] = 3
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.


シーケンスはリストから導き出すことができて，リストが使える文脈では常に利用することができる:

::

    sage: v = Sequence([1,2,3,4/5])
    sage: isinstance(v, list)
    True
    sage: list(v)
    [1, 2, 3, 4/5]
    sage: type(list(v))
    <... 'list'>


不変性シーケンスの例としては，ベクトル空間の基底系があげられる．
基底系そのものが変わっては困るから，これは当然のことだ．


::

    sage: V = QQ^3; B = V.basis(); B
    [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
    ]
    sage: type(B)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: B[0] = B[1]
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.
    sage: B.universe()
    Vector space of dimension 3 over Rational Field


ディクショナリ
===============

ディクショナリ(「連想配列」と呼ばれる場合もある)とは，文字列、数値、タプルなどのハッシュ可能なオブジェクトから任意のオブジェクトへの写像のことである．
(ハッシュ可能オブジェクトについての詳細は http://docs.python.org/tut/node7.html と http://docs.python.org/lib/typesmapping.html を参照．)

::

    sage: d = {1:5, 'sage':17, ZZ:GF(7)}
    sage: type(d)
    <... 'dict'>
    sage: list(d.keys())
    [1, 'sage', Integer Ring]
    sage: d['sage']
    17
    sage: d[ZZ]
    Finite Field of size 7
    sage: d[1]
    5

三番目の例を見ると分るように，ディクショナリのインデックス(キー)として整数環のような複雑なオブジェクトでも使うことができる．


上の例のディクショナリは，同じデータを含むリストに直すことができる:


.. link

::

    sage: list(d.items())
    [(1, 5), ('sage', 17), (Integer Ring, Finite Field of size 7)]

ディクショナリに含まれるキーと値の対を反復に利用する場合に，よく使われるイディオムがある:


::

    sage: d = {2:4, 3:9, 4:16}
    sage: [a*b for a, b in d.items()]
    [8, 27, 64]

最後の出力を見ると判るように，ディクショナリ内は整列されていない．


集合
=====

Pythonには集合(set)型が組込まれている．
集合型の主な利点としては，標準的な集合演算が可能になるだけではなく，ある要素が集合に属するかどうかを極めて高速に判定する機能を備えている点があげられる．

::

    sage: X = set([1,19,'a']);   Y = set([1,1,1, 2/3])
    sage: X  # random sort order
    {1, 19, 'a'}
    sage: X == set(['a', 1, 1, 19])
    True
    sage: Y
    {2/3, 1}
    sage: 'a' in X
    True
    sage: 'a' in Y
    False
    sage: X.intersection(Y)
    {1}

さらに，Sageは(Pythonの組み込み集合型を使って実装されたものも含まれる)独自の集合型を備えており，こちらにはSageに固有の付加機能がいくつか加えられている．
このSage独自の集合型を生成するには， ``Set(...)`` を使う．
例えば

::

    sage: X = Set([1,19,'a']);   Y = Set([1,1,1, 2/3])
    sage: X # random sort order
    {'a', 1, 19}
    sage: X == Set(['a', 1, 1, 19])
    True
    sage: Y
    {1, 2/3}
    sage: X.intersection(Y)
    {1}
    sage: print(latex(Y))
    \left\{1, \frac{2}{3}\right\}
    sage: Set(ZZ)
    Set of elements of Integer Ring


イテレータ
===========

イテレータ(iterator)は最近になってPythonに加えられた機能で，数学指向のアプリケーション作成にはとりわけ便利なものだ．
以下で実例を見ていくことにするが，使用法の詳細は [PyT]_ を見てほしい．
まず :math:`10000000` までの非負整数の平方に関するイテレータを作ってみよう．

::

    sage: v = (n^2 for n in range(10000000))
    sage: next(v)
    0
    sage: next(v)
    1
    sage: next(v)
    4

今度は，素数 :math:`p` から :math:`4p+1` の形の素数に関するイテレータを作り，最初の数個を見てみることにする．


::

    sage: w = (4*p + 1 for p in Primes() if is_prime(4*p+1))
    sage: w         # 次の行の 0xb0853d6c はランダムに生成された16進数
    <generator object <genexpr> at ...>
    sage: next(w)
    13
    sage: next(w)
    29
    sage: next(w)
    53

有限体，整数など，ある種の環にはイテレータが付随している:

::

    sage: [x for x in GF(7)]
    [0, 1, 2, 3, 4, 5, 6]
    sage: W = ((x,y) for x in ZZ for y in ZZ)
    sage: next(W)
    (0, 0)
    sage: next(W)
    (0, 1)
    sage: next(W)
    (0, -1)


ループ，関数，制御文，比較
=============================

``for`` ループの一般的な使用法については，これまでに何度も実例を見ている．
Pythonでは、 ``for`` ループ構文はインデントで分節されている．
次のような具合だ:

::

    >>> for i in range(5):
    ...     print(i)
    ...
    0
    1
    2
    3
    4


``for`` 文はコロン ``:`` で終っており，ループの本体すなわち ``print(i)`` がインデントされていることに注意(GAPやMapleに見られる "do"や "od"はない)．
Pythonでは，このインデントが重要な役割を果たしている．
以下の例のように，Sageでは ``:`` に続けて ``enter`` キーを押すと自動的にインデントが挿入される．

::

    sage: for i in range(5):
    ....:     print(i)  # ここでは[Enter]を2回押す
    ....:
    0
    1
    2
    3
    4


代入には ``=`` 記号，比較には ``==`` 記号を使う:

::

    sage: for i in range(15):
    ....:     if gcd(i,15) == 1:
    ....:         print(i)
    ....:
    1
    2
    4
    7
    8
    11
    13
    14

``if`` ， ``for`` および ``while`` 文のブロック構造が，インデントによって決まっているところに注目：


::

    sage: def legendre(a,p):
    ....:     is_sqr_modp=-1
    ....:     for i in range(p):
    ....:         if a % p == i^2 % p:
    ....:             is_sqr_modp=1
    ....:     return is_sqr_modp

    sage: legendre(2,7)
    1
    sage: legendre(3,7)
    -1

むろん，上のコードの目的はPython/Sageによるプログラムの特徴を例示することであって，Legendre記号の効率的な実装にはなっていない．
Sageに付属している関数 ``kronecker`` は，PARIのCライブラリを経由してLegendre記号を効率良く計算する．


最後に注意したいのは ``==`` ， ``!=`` ， ``<=`` ， ``>=`` ， ``>`` ， ``<`` などを使った比較演算では，比較対象の量は可能ならば自動的に同じ型に変換されることだ:


::

    sage: 2 < 3.1; 3.1 <= 1
    True
    False
    sage: 2/3 < 3/2;   3/2 < 3/1
    True
    True

記号を含む不等号の判定には  ``bool`` 関数を用いる:

::

    sage: x < x + 1
    x < x + 1
    sage: bool(x < x + 1)
    True


Sageにおける異種オブジェクト間の比較演算では，まず対象オブジェクトの共通ペアレント型への正準型強制(変換)が試みられる(「型強制」(coercion)の詳細については :ref:`section-coercion` 節を参照)．
比較演算は，この型強制が成功後，変換されたオブジェクトに対して実行される．
変換が不可能であれば，その二つのオブジェクトは等しくないと判定されることになる．
二つの変数が同一のオブジェクトを参照(レファレンス)しているかどうかを調べるには、 ``is`` を使う．
例えばPythonのint型 ``1`` は唯一だが，SageのInteger型 ``1`` は違う:

::

    sage: 1 is 2/2
    False
    sage: 1 is 1
    False
    sage: 1 == 2/2
    True

以下に示す例の，前半の二つの不等式は  ``False`` になる．
これは正準写像 :math:`\QQ\to \GF{5}` が存在せず， :math:`\GF{5}` 上の :math:`1` を :math:`1 \in \QQ` と比較する基準がないためである．
一方、後半の不等式は関係する正準写像 :math:`\ZZ \to \GF{5}` が存在するため ``True`` と判定される．
比較する式の左辺右辺を入れ替えても結果は変わらない．


::

    sage: GF(5)(1) == QQ(1); QQ(1) == GF(5)(1)
    False
    False
    sage: GF(5)(1) == ZZ(1); ZZ(1) == GF(5)(1)
    True
    True
    sage: ZZ(1) == QQ(1)
    True

*警告:* Sageにおける比較演算は， :math:`1 \in \GF{5}` は :math:`1 \in  \QQ` と等しいとみなすMagmaよりも制限がきつい．

::

    sage: magma('GF(5)!1 eq Rationals()!1')  # optional - magma オプションでmagmaが必要
    true


プロファイリング
================

著者: Martin Albrecht (malb@informatik.uni-bremen.de)

    「早計な最適化は，あらゆる災厄の源である」 -- ドナルド・クヌース


コードのボトルネック，つまり処理時間の大半を費している部分を洗い出さなければならない場面がある．
ボトルネックが判らなければどこを最適化すべきかの判定もできないからだ．
コードのボトルネックを特定する作業のことをプロファイリングと呼ぶが，Python/Sageにはプロファイリングの手段が何通りも用意されている．


一番簡単な方法は、対話型シェルで ``prun`` コマンドを実行することだ．
``prun``  は，使われた関数と各関数の処理時間の一覧を表示してくれる．
例として，有限体上の行列積演算(バージョン1.0で、まだ低速だ)をプロファイルしてみよう．その手順は:


::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])

.. skip


::

    sage: %prun B = A*A
           32893 function calls in 1.100 CPU seconds

    Ordered by: internal time

    ncalls tottime percall cumtime percall filename:lineno(function)
     12127  0.160   0.000   0.160  0.000 :0(isinstance)
      2000  0.150   0.000   0.280  0.000 matrix.py:2235(__getitem__)
      1000  0.120   0.000   0.370  0.000 finite_field_element.py:392(__mul__)
      1903  0.120   0.000   0.200  0.000 finite_field_element.py:47(__init__)
      1900  0.090   0.000   0.220  0.000 finite_field_element.py:376(__compat)
       900  0.080   0.000   0.260  0.000 finite_field_element.py:380(__add__)
         1  0.070   0.070   1.100  1.100 matrix.py:864(__mul__)
      2105  0.070   0.000   0.070  0.000 matrix.py:282(ncols)
      ...

ここで ``ncalls`` は関数の呼出し回数， ``tottime`` はその関数が費した総時間(配下の関数群の呼出しにかかった時間は除く)， ``percall`` は ``tottime`` を ``ncalls`` で割って得られる平均総消費時間， ``cumtime``  は関数本体とそれが呼出している配下の全関数双方の処理にかかった(つまり注目している関数の呼出しから開放までの)全時間， ``percall`` は ``cumtime`` を呼出し回数で割って得られる関数の平均処理時間で， ``filename:lineno(function)`` は各関数の所在情報を示している．
結局のところ， ``prun`` の表示一覧の上位にある関数ほど処理に要する負荷も大きいことが判る．
これで，どこを最適化すべきか考えやすくなるはずだ．


これまでと同じように， ``prun?`` と入力するとプロファイラの使い方と表示情報の解釈について詳細を見ることができる．


プロファイル情報をオブジェクトとして保存しておいて，後で詳しく調べてもよい:


.. skip

::

    sage: %prun -r A*A
    sage: stats = _
    sage: stats?

*注意*: ``stats = prun -r A\*A`` と実行すると，文法エラーが表示される．
これは ``prun`` がIPythonのシェルコマンドであって，通常の関数ではないためである．


プロファイル情報をグラフィカルに表示したければ，hotshotプロファイラ， ``hotshot2cachetree`` スクリプトと ``kcachegrind`` プログラム(Unix系のみ)などを使えばよい．
上の例と同じ作業をhotshotプロファイラで実行すると:

.. skip

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])
    sage: import hotshot
    sage: filename = "pythongrind.prof"
    sage: prof = hotshot.Profile(filename, lineevents=1)

.. skip

::

    sage: prof.run("A*A")
    <hotshot.Profile instance at 0x414c11ec>
    sage: prof.close()

得られた結果は，現ディレクトリのファイル ``pythongrind.prof`` に保存されている．
これをcachegrindフォーマットに変換すれば，内容をビジュアル化して把握することができる．


システムのコマンドシェルに戻り

.. skip

::

    hotshot2calltree -o cachegrind.out.42 pythongrind.prof

として変換を実行すると，出力ファイル ``cachegrind.out.42`` 内の情報を ``kcachegrind`` コマンドを使って検討することができる．
ファイルの慣例的名称 ``cachegrind.out.XX`` は残念ながら変更できない．

