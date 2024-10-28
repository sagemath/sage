.. _chapter-help:

ヘルプの利用
============

Sageには充実したドキュメントが組込まれていて，(例えば)関数や定数の名前に続けて疑問符 ``?`` を入力するだけで解説を呼び出すことができる:

.. skip

::

    sage: tan?
    Type:        <class 'sage.calculus.calculus.Function_tan'>
    Definition:  tan( [noargspec] )
    Docstring:

        The tangent function

        EXAMPLES:
            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
    sage: log2?
    Type:        <class 'sage.functions.constants.Log2'>
    Definition:  log2( [noargspec] )
    Docstring:

        The natural logarithm of the real number 2.

        EXAMPLES:
            sage: log2
            log2
            sage: float(log2)
            0.69314718055994529
            sage: RR(log2)
            0.693147180559945
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(log2)
            0.69314718055994530941723212145817656807550013436025525412068
            sage: l = (1-log2)/(1+log2); l
            (1 - log(2))/(log(2) + 1)
            sage: R(l)
            0.18123221829928249948761381864650311423330609774776013488056
            sage: maxima(log2)
            log(2)
            sage: maxima(log2).float()
            .6931471805599453
            sage: gp(log2)
            0.6931471805599453094172321215             # 精度は32ビット
            0.69314718055994530941723212145817656807   # 同じく64ビット
    sage: sudoku?
    File:        sage/local/lib/python2.5/site-packages/sage/games/sudoku.py
    Type:        <... 'function'>
    Definition:  sudoku(A)
    Docstring:

        Solve the 9x9 Sudoku puzzle defined by the matrix A.

        EXAMPLE:
            sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0,
        0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0, 0,0,0, 2,0,8, 0,0,0,
        0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2,
        0,0,0, 4,9,0, 0,5,0, 0,0,3])
            sage: A
            [5 0 0 0 8 0 0 4 9]
            [0 0 0 5 0 0 0 3 0]
            [0 6 7 3 0 0 0 0 1]
            [1 5 0 0 0 0 0 0 0]
            [0 0 0 2 0 8 0 0 0]
            [0 0 0 0 0 0 0 1 8]
            [7 0 0 0 0 4 1 5 0]
            [0 3 0 0 0 2 0 0 0]
            [4 9 0 0 5 0 0 0 3]
            sage: sudoku(A)
            [5 1 3 6 8 7 2 4 9]
            [8 4 9 5 2 1 6 3 7]
            [2 6 7 3 4 9 5 8 1]
            [1 5 8 4 6 3 9 7 2]
            [9 7 4 2 1 8 3 6 5]
            [3 2 6 7 9 5 4 1 8]
            [7 8 2 9 3 4 1 5 6]
            [6 3 5 1 7 2 8 9 4]
            [4 9 1 8 5 6 7 2 3]

さらにSageでは，関数名の最初の何文字かを入力して ``TAB`` キーを打つと候補が表示される，タブ補完機能を使うことができる．
例えば ``ta`` と入力して ``TAB`` キーを押せば、Sageは ``tachyon, tan, tanh,taylor`` を表示してくる．
この機能を使えば関数やオブジェクトなどの名称を容易に探しあてることができるはずだ．


.. _section-functions:

関数, インデントおよび数え上げ
====================================

Sageで新しい関数を定義するには， ``def`` 命令を使い、変数名を並べた後にコロン ``:`` を付ける．
以下に例を示そう:

::

    sage: def is_even(n):
    ....:     return n % 2 == 0
    ....:
    sage: is_even(2)
    True
    sage: is_even(3)
    False

*注意* : チュートリアルをどの形式で閲覧しているかにもよるが，上のコード例の2行目には四つのドット ``....`` が見えているはずだ．
この四点ドットは入力しないこと．
四点ドットは，コードがインデントされていることを示しているだけだからだ．
そうした場面では，常に構文ブロックの末尾で一度 ``Return/Enter`` を押して空行を挿入し，関数定義を終了してやらねばならない．

引数の型を指定していないことに注意．複数個の引数を指定し，その各々にデフォルト値を割り当てることもできる．
例えば、以下の関数では引数 ``divisor`` の値が指定されない場合， ``divisor=2``  がデフォルト値になる:

::

    sage: def is_divisible_by(number, divisor=2):
    ....:     return number % divisor == 0
    sage: is_divisible_by(6,2)
    True
    sage: is_divisible_by(6)
    True
    sage: is_divisible_by(6, 5)
    False


関数を呼び出すときには，特定の引数へ明示的に値を代入することもできる．
引数への明示的な代入を行なう場合，関数に渡す引数の順序は任意になる:

.. link

::

    sage: is_divisible_by(6, divisor=5)
    False
    sage: is_divisible_by(divisor=2, number=6)
    True


Pythonの構文ブロックは，他の多くの言語のように中括弧やbegin-endで括ることによって示されるわけではない．
代りに、Pythonでは構文構造に正確に対応したインデンテーション(字下げ)によってブロックを示す．
次の例は， ``return`` ステートメントが関数内の他のコードと同じようにインデントされていないために，文法エラーになっている:

.. skip

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3, n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:    return v
    Syntax Error:
           return v

しかし正しくインデントし直せば，この関数はきちんと動くようになる:

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:     return v
    sage: even(10)
    [4, 6, 8]

行末にセミコロンは必要ない．
ほとんどの場合，行末は改行記号によって示される．
しかし，1行に複数のステートメントをセミコロンで区切って書き込むこともできる:

::

    sage: a = 5; b = a + 3; c = b^2; c
    64

1行の内容を複数行に分けて書きたければ，各行末にバックスラッシュをつければよい:

::

    sage: (2 +
    ....:    3)
    5


Sageでは，一定範囲の整数の数え上げによって反復を制御する．
例えば，以下のコードの1行目はC++やJavaにおける ``for(i=0; i<3; i++)``  と全く同じ意味になる:


::

    sage: for i in range(3):
    ....:     print(i)
    0
    1
    2


次の例の最初の行は， ``for(i=2;i<5;i++)`` に対応している．


::

    sage: for i in range(2,5):
    ....:     print(i)
    2
    3
    4


``range`` の三つ目の引数は増分値を与えるので,
次のコードは ``for(i=1;i<6;i+=2)`` と同じ意味になる.

::

    sage: for i in range(1,6,2):
    ....:     print(i)
    1
    3
    5


Sageで計算した値を見映えよく表形式に並べて表示したくなることもあるだろう．
そんなとき役立つのが文字列フォーマットだ．
以下では，各列幅がきっかり6文字分の表を作り，整数とその2乗、3乗の値を並べてみる．


::

    sage: for i in range(5):
    ....:     print('%6s %6s %6s' % (i, i^2, i^3))
         0      0      0
         1      1      1
         2      4      8
         3      9     27
         4     16     64


Sageにおける最も基本的なデータ構造はリストで，名前の示すとおり任意のオブジェクトの並びのことである．
上で使った ``range`` も、整数のリストを生成している：

::

    sage: list(range(2,10))
    [2, 3, 4, 5, 6, 7, 8, 9]

もう少し複雑なリストの例として:

::

    sage: v = [1, "hello", 2/3, sin(x^3)]
    sage: v
    [1, 'hello', 2/3, sin(x^3)]

多くのプログラミング言語と同じように，リスト添字は0から始まる．


.. link

::

    sage: v[0]
    1
    sage: v[3]
    sin(x^3)


``v`` の長さを取得するには ``len(v)`` を使い， ``v`` の末尾に新しいオブジェクトを追加するには ``v.append(obj)`` ，そして ``v`` の :math:`i` 番目の要素を削除するためには ``del v[i]`` とする:


.. link


::

    sage: len(v)
    4
    sage: v.append(1.5)
    sage: v
    [1, 'hello', 2/3, sin(x^3), 1.50000000000000]
    sage: del v[1]
    sage: v
    [1, 2/3, sin(x^3), 1.50000000000000]


もう一つの重要なデータ構造がディクショナリ(連想配列とも言う)である．
ディクショナリの振舞いはリストに似ているが，異なるのはその添字付けに基本的にいかなるオブジェクトでも使うことができる点だ(ただし添字は不変性オブジェクトでなければならない)．

::

    sage: d = {'hi':-2,  3/8:pi,   e:pi}
    sage: d['hi']
    -2
    sage: d[e]
    pi


クラスを使えば自分で新しいデータ型を定義することも可能だ．
クラスによる数学オブジェクトのカプセル化は，Sageプログラムを見通しよく構成するための強力な方法である．
以下では， *n* までの正の偶数のリストを表すクラスを定義してみよう．
定義には組み込み型 ``list`` を使っている．
::

    sage: class Evens(list):
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:         list.__init__(self, range(2, n+1, 2))
    ....:     def __repr__(self):
    ....:         return "Even positive numbers up to n."


オブジェクトの生成時には初期化のために ``__init__`` メソッドが呼ばれ， ``__repr__`` メソッドはオブジェクトを印字する．
``__init__`` メソッドの2行目ではリストコンストラクタを使った．
クラス ``Evens`` のオブジェクトを生成するには以下のようにする:


.. link


::

    sage: e = Evens(10)
    sage: e
    Even positive numbers up to n.


``e`` の印字には，我々が定義した ``__repr__`` メソッドが使われている．
オブジェクトに含まれる偶数のリストを表示するには， ``list`` 関数を使う:


.. link


::

    sage: list(e)
    [2, 4, 6, 8, 10]


``n`` 属性にアクセスし、 ``e`` をリストのように扱うこともできる．


.. link


::

    sage: e.n
    10
    sage: e[2]
    6
