.. _chapter-interactive_shell:

*********************
対話型シェル
*********************

このチュートリアルの大部分は，読者が ``sage`` コマンドによってSageインタプリタを起動しているものと前提している．
コマンド ``sage`` は改造版IPythonシェルを起動し，大量の関数やクラス群をインポートしてコマンドプロンプトから利用可能にする．
``$SAGE_ROOT/ipythonrc`` ファイルを編集すれば，さらなるシェル環境のカスタマイズも可能だ．
Sageを起動すると，すぐに次のような画面が現れる:

.. skip

::

    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘


    sage:

Sageを終了するには，Ctrl-Dと押すか， コマンド ``quit`` あるいは ``exit`` を入力する．

.. skip

::

    sage: quit
    Exiting Sage (CPU time 0m0.00s, Wall time 0m0.89s)


"Wall time"は，CPUタイムではなく外界の実経過時間を示している．
CPUタイムはGAPやSingularなどのサブプロセスの消費時間までは勘定に入れてくれないから，実経過時間も計算時間の見積りに必要だ．

(ターミナルから ``kill -9`` を入力してSageプロセスを停止するのは止めたほうがいい．
``kill -9`` ではMapleなどの子プロセスが停止しなかったり， ``$HOME/.sage/tmp`` 内の一時ファイルが消去されずに終わるなどの恐れがある．)


Sageセッション
=================

*セッション* とは，Sageの起動から終了までの間に行なわれた一連の入出力の総体のことをいう．
Sageは，Sageに対する入力の全てをIPython経由で記録している．
事実，(ノートブック経由ではなく)対話型シェルを使ってSageを動かしているのならば，好きな時に ``%history`` (または ``%hist``) と入力して，それまでの全入力履歴を見ることができる．
IPythonについてもっと知りたければ、Sageプロンプトで ``?`` と入力すれば， "IPython offers numbered prompts ... with input and output caching. All input is saved and can be retrieved as variables (besides the usual arrow key recall). The following GLOBAL variables always exist (so don't overwrite them!)" などと詳しい情報を表示させることができる:

::

       _:  前回の入力を呼び出す (対話型シェルとノートブックの両方で通用する)
      __: 前々回の入力を呼び出す(対話型シェルのみで通用)
    _oh : 全ての入力をリストする(対話型シェルのみで通用)

ここで例を見てみよう:

.. skip

::

    sage: factor(100)
     _1 = 2^2 * 5^2
    sage: kronecker_symbol(3,5)
     _2 = -1
    sage: %hist   # これが使えるのは対話型シェル上のみ．ノートブックではだめ．
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    sage: _oh
     _4 = {1: 2^2 * 5^2, 2: -1}
    sage: _i1
     _5 = 'factor(ZZ(100))\n'
    sage: eval(_i1)
     _6 = 2^2 * 5^2
    sage: %hist
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    4: _oh
    5: _i1
    6: eval(_i1)
    7: %hist

以降，このチュートリアル，それに他のSageドキュメンテーションでも出力番号を省略する．

セッション中は，一連の入力をマクロとして保存しておいて再利用することもできる．


.. skip

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: M = ModularSymbols(37)
    sage: %hist
    1: E = EllipticCurve([1,2,3,4,5])
    2: M = ModularSymbols(37)
    3: %hist
    sage: %macro em 1-2
    Macro `em` created. To execute, type its name (without quotes).


.. skip

::

    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E = 5
    sage: M = None
    sage: em
    Executing Macro...
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field

対話型シェルを使っている間も，感嘆符 ``!`` を前置すれば好きなUNIXシェルコマンドを実行することができる．
例えば

.. skip

::

    sage: !ls
    auto  example.sage glossary.tex  t  tmp  tut.log  tut.tex

のように，カレントディレクトリの内容を表示することができる．

シェル変数 ``PATH`` の先頭にはSageのbinディレクトリが配置されているから， ``gp`` ， ``gap`` ， ``singular`` ， ``maxima`` などを実行すると，Sageに付属しているプログラムのバージョンを確認することができる．

.. skip


::

    sage: !gp
    Reading GPRC: /etc/gprc ...Done.

                               GP/PARI CALCULATOR Version 2.2.11 (alpha)
                      i686 running linux (ix86/GMP-4.1.4 kernel) 32-bit version
    ...
    sage: !singular
                         SINGULAR                             /  Development
     A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \


入出力のログをとる
========================

Sageセッションのロギングと，セッションの保存(:ref:`section-save` 節を参照)は同じことではない．
入力のログをとるには， ``logstart`` コマンドを使う(オプションで出力のログも可能だ)．
詳細については ``logstart?`` と入力してみてほしい．
``logstart`` を使えば，全ての入力と出力のログを残し，将来のセッション時に(そのログファイルをリロードしてやるだけで)入力を再生することも可能になる．

.. skip

::

    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: logstart setup
    Activating auto-logging. Current session state plus future input saved.
    Filename       : setup
    Mode           : backup
    Output logging : False
    Timestamping   : False
    State          : active
    sage: E = EllipticCurve([1,2,3,4,5]).minimal_model()
    sage: F = QQ^3
    sage: x,y = QQ['x,y'].gens()
    sage: G = E.gens()
    sage:
    Exiting Sage (CPU time 0m0.61s, Wall time 0m50.39s).
    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: load("setup")
    Loading log file <setup> one line at a time...
    Finished replaying log file <setup>
    sage: E
    Elliptic Curve defined by y^2 + x*y  = x^3 - x^2 + 4*x + 3 over Rational
    Field
    sage: x*y
    x*y
    sage: G
    [(2 : 3 : 1)]

SageをLinux KDEターミナル ``konsole`` 上で使っているなら，以下の手順でセッションを保存することもできる．
まず ``konsole`` 上でSageを起動したら、 "settings"(日本語環境であれば『設定』)を選択し，次に "history"(『履歴』)， "set unlimited"(『無制限にする』)の順に選択しておく．
セッションを保存したくなった時点で， "edit"(『編集』)の中の "save history as..."(『履歴を名前を付けて保存』)を選択してセッションを保存するファイル名を入力してやればよい．
いったんファイルとして保存してしまえば，好きなようにxemacsなどのエディタで読み込んだりプリントアウトしたりすることができる．


プロンプト記号はペースト時に無視される
========================================

SageセッションあるいはPythonの演算結果を読み込んで，Sage上にコピーしたい場合がある．
厄介なのは、そうした出力に ``>>>`` や ``sage:`` といったプロンプト記号が紛れ込んでいることだ．
しかし実際には，プロンプト記号を含む実行例をSage上へ好きにコピー・ペーストしてやることができる．
デフォルトでSageパーサーはデータをPythonに送る前に行頭の ``>>>`` や ``sage:`` プロンプト記号を除去してくれるからだ．
例えば

.. skip

::

    sage: 2^10
    1024
    sage: sage: sage: 2^10
    1024
    sage: >>> 2^10
    1024


計時コマンド
===============

入力行の先頭に ``%time`` コマンドを入れておくと，出力までに要した時間を表示することができる．
例として，べき乗計算を異なった方法で行なった場合の実行時間を比較してみよう．
以下に示した実行時間の値は，動かしているコンピュータ本体やSageのバージョンによって大きく異なる可能性が高い．
まず、Pythonを直に動かしてみると:

.. skip

::

    sage: %time a = int(1938)^int(99484)
    CPU times: user 0.66 s, sys: 0.00 s, total: 0.66 s
    Wall time: 0.66


上の出力は，実行に計0.66秒かかり， "Wall time" つまりユーザーの実待ち時間もやはり0.66秒だったことを示している．
コンピュータに他のプログラムから大きな負荷がかかっている場合， "Wall time"がCPUタイムよりかなり長くなることがある．

次に，同じべき乗計算をSage組み込みのInteger型を使って実行した場合の時間を計ってみよう．
SageのInteger型は，Cython経由でGMPライブラリを使って実装されている:

.. skip

::

    sage: %time a = 1938^99484
    CPU times: user 0.04 s, sys: 0.00 s, total: 0.04 s
    Wall time: 0.04

PARIのC-ライブラリを経由すると

.. skip

::

    sage: %time a = pari(1938)^pari(99484)
    CPU times: user 0.05 s, sys: 0.00 s, total: 0.05 s
    Wall time: 0.05

GMPの方が速いが，その差はわずかだ(Sage用にビルドされたPARIは整数演算にGMPを使っているのだから，納得できる結果である)．


次の例のように， ``cputime`` コマンドを使えば，一連のコマンドからなるコードブロックの実行時間を計ることもできる:

::

    sage: t = cputime()
    sage: a = int(1938)^int(99484)
    sage: b = 1938^99484
    sage: c = pari(1938)^pari(99484)
    sage: cputime(t)                       # random 値には若干の幅がある．
    0.64


.. skip


::

   sage: cputime?
    ...
        Return the time in CPU second since Sage started, or with optional
        argument t, return the time since time t.
        INPUT:
            t -- (optional) float, time in CPU seconds
        OUTPUT:
            float -- time in CPU seconds


``walltime`` コマンドの動作は，計測するのが実経過時間である点以外は ``cputime`` コマンドと変わらない．

上で求めたべき乗を，Sageに取り込まれている各コンピュータ代数システムを使って計算することもできる．
計算を実行するには，使いたいシステムの名前をコマンド名としてそのプログラムのサーバを呼び出す．
いちばん肝心な計測値は，実経過時間(wall time)だ．
しかし，実経過時間とCPUタイムの値が大幅に食い違う場合は，解決すべきパフォーマンス上の問題点の存在を示している可能性がある．

.. skip

::

    sage: time 1938^99484;
    CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
    Wall time: 0.01
    sage: gp(0)
    0
    sage: time g = gp('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: maxima(0)
    0
    sage: time g = maxima('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.30
    sage: kash(0)
    0
    sage: time g = kash('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: mathematica(0)
            0
    sage: time g = mathematica('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.03
    sage: maple(0)
    0
    sage: time g = maple('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.11
    sage: gap(0)
    0
    sage: time g = libgap.eval('1938^99484;')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 1.02

以上のテスト計算で最も遅かったのは，GAPとMaximaである(実行結果はホスト ``sage.math.washington.edu`` 上のもの)．
各システムとのpexpectインターフェイスにかかる負荷を考えると，上の一連の計測値を最速だったSageの値と比較するのは公平を欠く面があるかもしれない．


IPythonトリック
====================

すでに述べたように，SageはそのフロントエンドとしてIPythonを援用しており，ユーザはIPythonのコマンドと独自機能を自由に利用することができる．
その全貌については， ご自分で `full IPython documentation
<http://ipython.scipy.org/moin/Documentation>`_ を読んみてほしい．
そのかわり，ここではIPythonの「マジックコマンド」と呼ばれる，お便利なトリックをいくつか紹介させていただこう:

- ``%edit`` (``%ed`` や ``ed`` でもいい)を使ってエディタを起動すれば，複雑なコードの入力が楽になる．
  Sageの使用前に，環境変数 :envvar:`EDITOR` に好みのエディタ名を設定しておこう(``export EDITOR=/usr/bin/emacs`` または ``export EDITOR=/usr/bin/vim`` とするか， ``.profile`` ファイルなどで同様の設定をする)．
  するとSageプロンプトで ``%edit`` を実行すれば設定したエディタが起動する．そのエディタで関数

  ::

    def some_function(n):
        return n**2 + 3*n + 2

  を定義し，保存したらエディタを終了する．
  以降，このセッション中は ``some_function`` を利用できるようになる．
  内容を編集したければSageプロンプトで ``%edit some_function`` と入力すればよい．


- 結果出力を他の用途のために編集したければ， ``%rep`` を実行する．
  すると直前に実行したコマンドの出力が編集できるようにSageプロンプト上に配置される．::

   sage: f(x) = cos(x)
   sage: f(x).derivative(x)
   -sin(x)

  この段階でSageプロンプトから ``%rep`` を実行すると，新しいSageプロンプトに続いて ``-sin(x)`` が現われる．
  カーソルは同じ行末にある．


IPythonのクイック レファレンスガイドを見たければ， ``%quickref`` と入力する．
執筆時点(2011年4月)ではSageはIPythonのバージョン0.9.1を採用しており， `documentation for its magic commands
<http://ipython.org/ipython-doc/dev/interactive/tutorial.html#magic-functions>`_
はオンラインで読むことができる．
マジックコマンドの，ちょっと進んだ機能群についてはIPythonの `ここ
<http://ipython.org/ipython-doc/stable/interactive/reference.html#magic-command-system>`_
で文書化されているのが見つかるはずだ．


エラーと例外処理
=====================

処理中に何かまずいことが起きると，Pythonはふつう『例外』(exception)を発生し，その例外を引き起こした原因を教えてくれることもある．
よくお目にかかることになるのは， :class:`NameError` や :class:`ValueError` といった名称の例外だ(Pythonライブラリーリファレンス [PyLR]_ に例外名の包括的なリストがある)．
実例を見てみよう:

::

    sage: EllipticCurve([0,infinity])
    Traceback (most recent call last):
    ...
    SignError: cannot multiply infinity by zero


何が悪いか調べるには対話型デバッガが役立つこともある．
デバッガを使うには、 ``%pdb`` コマンドによって作動のオン/オフをトグルする(デフォルトはオフ)．
作動後は、例外が発生するとデバッガが起動し，プロンプト ``ipdb>`` が表示される．
このデバッガの中から，任意のローカル変数の状態を表示したり，実行スタックを上下して様子を調べることができる．

.. skip

::

    sage: %pdb
    Automatic pdb calling has been turned ON
    sage: EllipticCurve([1,infinity])
    ---------------------------------------------------------------------------
    <class 'exceptions.TypeError'>             Traceback (most recent call last)
    ...

    ipdb>


デバッガから実行できるコマンドの一覧を見るには， ``ipdb>`` プロンプト上で ``?`` を入力する:
::

    ipdb> ?

    Documented commands (type help <topic>):
    ========================================
    EOF    break  commands   debug    h       l     pdef   quit    tbreak
    a      bt     condition  disable  help    list  pdoc   r       u
    alias  c      cont       down     ignore  n     pinfo  return  unalias
    args   cl     continue   enable   j       next  pp     s       up
    b      clear  d          exit     jump    p     q      step    w
    whatis where

    Miscellaneous help topics:
    ==========================
    exec  pdb

    Undocumented commands:
    ======================
    retval  rv

Sageに戻るには，Ctrl-Dか ``quit`` を入力する．


.. _section-tabcompletion:

コマンド入力の遡行検索とタブ補完
=================================

*遡行検索*: コマンドの冒頭部を打ち込んでから ``Ctrl-p``  (または上向き矢印キー)を押すと，冒頭部が一致する過去の入力行を全て呼び出すことができる．
この機能は，Sageをいったん終了し再起動してからでも有効である．
``Ctrl-r`` を入力すれば，入力ヒストリを逆方向に検索することも可能だ．
この入力行の検索と再利用機能は全て ``readline``  パッケージを経由しており，ほとんどのLinux系システム上で利用できるはずだ．

タブ補完機能を体験するため，まず3次元ベクトル空間 :math:`V=\QQ^3` を生成しておく:
::

    sage: V = VectorSpace(QQ,3)
    sage: V
    Vector space of dimension 3 over Rational Field

次のような，もっと簡潔な記号法を使ってもよい:

::

    sage: V = QQ^3


タブ補完を使えば，簡単に :math:`V` の全メンバ関数を一覧表示することができる．
``V.`` と入力し，ついで ``[tab]`` キーを押すだけだ:

.. skip

::

    sage: V.[tab key]
    V._VectorSpace_generic__base_field
    ...
    V.ambient_space
    V.base_field
    V.base_ring
    V.basis
    V.coordinates
    ...
    V.zero_vector

関数名の出だし何文字かを打ってから ``[tab キー]`` を押せば，入力した文字で始まる名前の関数だけに候補を絞ることができる．

.. skip

::

    sage: V.i[tab key]
    V.is_ambient  V.is_dense    V.is_full     V.is_sparse

特定の関数について調べたい場合もある．
coordinates関数を例にとると，そのヘルプを表示するには ``V.coordinates?`` と入力すればいいし，ソースコードを見るには ``V.coordinates??`` を入力すればいい．
詳細については次の節で解説する．


統合ヘルプシステム
======================

Sageの特長の一つは，総合的なヘルプ機能の装備である．
関数名に続けて?を入力すると、その関数のドキュメントを表示することができる．

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates?
    Type:           instancemethod
    Base Class:     <class 'instancemethod'>
    String Form:    <bound method FreeModule_ambient_field.coordinates of Vector
    space of dimension 3 over Rational Field>
    Namespace:      Interactive
    File:           /home/was/s/local/lib/python2.4/site-packages/sage/modules/f
    ree_module.py
    Definition:     V.coordinates(self, v)
    Docstring:
        Write v in terms of the basis for self.

        Returns a list c such that if B is the basis for self, then

                sum c_i B_i = v.

        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: M = FreeModule(IntegerRing(), 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0-M1)
            [2, -1]


上で見たように，ヘルプ表示には，そのオブジェクトの型，定義されているファイル，現セッションにペーストすることができる使用例付きの解説が含まれる．
使用例のほとんどは常に自動的なテストが行なわれていて，仕様どおりの正確な動作が確認されている．

もう一つの機能は，Sageのオープンソース精神をよく表すものだ．
``f`` がPythonで書かれた関数であれば ``f??`` と入力すると ``f`` を定義しているソースを表示することができるのだ．
例えば

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates??
    Type:           instancemethod
    ...
    Source:
    def coordinates(self, v):
            """
            Write $v$ in terms of the basis for self.
            ...
            """
            return self.coordinate_vector(v).list()

これを見ると， ``coordinates`` 関数は ``coordinate_vector`` 関数を呼び出して結果をリストに変換しているだけであることが判る．
では ``coordinate_vector`` 関数が何をしているかと言うと:

.. skip

::

    sage: V = QQ^3
    sage: V.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            ...
            return self.ambient_vector_space()(v)


``coordinate_vector`` 関数は，入力を生成空間(ambient space)に合わせて型変換するから，これは :math:`v` の係数ベクトルが空間 :math:`V` ではどう変換されるか計算していることと同じである．
:math:`V` は :math:`\QQ^3` そのものだから，すでに同じ構造になっている．
部分空間用に，上とは異なる ``coordinate_vector`` 関数も用意されている．
部分空間を作って，どんな関数か見てみることにしよう:

.. skip

::

    sage: V = QQ^3; W = V.span_of_basis([V.0, V.1])
    sage: W.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            """
             ...
            """
            # First find the coordinates of v wrt echelon basis.
            w = self.echelon_coordinate_vector(v)
            # Next use transformation matrix from echelon basis to
            # user basis.
            T = self.echelon_to_user_matrix()
            return T.linear_combination_of_rows(w)

(こうした実装の仕方は無駄が多いと思われる方は，どうか我々に連絡して線形代数周りの最適化に力を貸していただきたい．)


``help(コマンド名)`` あるいは ``help(クラス名)`` と入力すれば，知りた いクラスのmanページ型ヘルプファイルを表示することもできる．


.. skip

::

    sage: help(VectorSpace)
    Help on function VectorSpace in module sage.modules.free_module:

    VectorSpace(K, dimension_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
                with_basis='standard', dimension=None, basis_keys=None, **args)
    EXAMPLES:

    The base can be complicated, as long as it is a field.

    ::

        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x
         over Integer Ring
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
    --More--


``q`` と入力してヘルプを終えると，中断前のセッション画面がそのまま復帰する．
セッションに干渉することがある ``function_name?`` と違って， ヘルプ表示はセッションの邪魔をしない．
とりわけ便利なのは  ``help(モジュール名)`` と入力することだ．
例えばベクトル空間は  ``sage.modules.free_module`` で定義されているから，そのモジュール全体に関するドキュメントを見たければ ``help(sage.modules.free_module)`` と実行すればよい．
ヘルプを使ってドキュメントを閲覧している間は， ``/`` と打てば語句検索 ができるし， ``?`` と打てば逆方向に検索することができる．


オブジェクトの保存と読み込み
=====================================

行列や，あるいはもっと手間のかかる複雑なモジュラーシンボルの空間を扱っていて，後で利用するため結果を保存しておきたくなったとしよう．
そんな場合にはどうすればよいだろうか．
オブジェクトを保存するために各コンピュータ代数システムが提供している方法は，以下の通りである．


#. **セッションの保存:** セッション全体の保存と読み込みのみ可能(GAP，Magmaなど)．

#. **統合入出力:** 全オブジェクトの印字が再利用可能な形式で行なわれる(GAPとPARI)．

#. **再実行**: インタープリタによるプログラムの再実行が容易にしてある(Singular，PARI)．


..
   #. **Save your Game:** Only support saving and loading of complete
      sessions (e.g., GAP, Magma).

   #. **Unified Input/Output:** Make every object print in a way that
      can be read back in (GP/PARI).

   #. **Eval**: Make it easy to evaluate arbitrary code in the
      interpreter (e.g., Singular, PARI).

Pythonで動くSageでは，全てのオブジェクトのシリアル化(直列化)という，他とは異なる方法が採用されている．
つまりオブジェクトを，その原型を再現可能な形式で文字列に変換するのだ．
これはPARIの統合入出力の考え方に近いが，オブジェクトを複雑な印字形式で画面出力してやる必要がないのが利点だ．
さらに保存と読み込みは(ほとんどの場合)完全に自動化されているから，新たにプログラムを書く必要もない．
そうした機能はPythonに最初から組込まれているものだからである．


ほぼ全てのSageオブジェクト ``x`` は， コマンド ``save(x,ファイル名)`` (あるいは多くの場合 ``x.save(ファイル名)``)を使えば圧縮形式でディスクに保存することができるようになっている．
保存したオブジェクトを読み戻すには， ``load(ファイル名)`` を実行する．


.. skip

::

    sage: A = MatrixSpace(QQ,3)(range(9))^2
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]
    sage: save(A, 'A')

ここでいったんSageを終了してみよう．再起動後に ``A``  を読み込むには:


.. skip

::

    sage: A = load('A')
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]

楕円曲線のようなもっと複雑なオブジェクトに対しても，以上と同じやり方が通用する．
メモリ上に配置されていたオブジェクト関連の全データは，そのオブジェクトと共に保存される．
例えば

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: v = E.anlist(100000)              # ちょっと時間がかかる
    sage: save(E, 'E')
    sage: quit

こうして保存された ``E`` は，オブジェクト本体と一緒に :math:`a_n` の冒頭100000個も保存するため，153Kバイトの大きさになる．


.. skip

::

    ~/tmp$ ls -l E.sobj
    -rw-r--r--  1 was was 153500 2006-01-28 19:23 E.sobj
    ~/tmp$ sage [...]
    sage: E = load('E')
    sage: v = E.anlist(100000)              # すぐ終了

(Python経由の保存と読み込みには， ``cPickle`` モジュールが使われている．
実際，Sageオブジェクト ``x`` の保存は ``cPickle.dumps(x, 2)`` を実行して行なうことができる．引数 ``2`` に注目．)


Sageで保存・読み込みできないのは，GAP， Singular， Maximaなど外部コンピュータ代数システムで作成されたオブジェクトである．
これらは読み込むことができても "invalid"(利用不能)な状態にあると認識される．
GAPでは，相当数のオブジェクトが再構成に使える印字形式を持つ一方，再構成できな場合も多いため印字形式からのオブジェクトの再構成は意図的に禁止されている．


.. skip

::

    sage: a = gap(2)
    sage: a.save('a')
    sage: load('a')
    Traceback (most recent call last):
    ...
    ValueError: The session in which this object was defined is no longer
    running.

GP/PARIオブジェクトは，印字形式から十分に再構成可能なため，保存と読み込みも可能になっている．


.. skip

::

    sage: a = gp(2)
    sage: a.save('a')
    sage: load('a')
    2


保存したオブジェクトは，異なるアーキテクチャ上の，異なるオペレーティングシステムで動くSageへもロードすることができる．
例えば，32ビット OSX上で保存した大規模行列を64ビット版LinuxのSageへ読み込み，その階段形式を求めてから元のOS X上へ戻すといったことも可能だ．
さらに，オブジェクトの保存までに使ったのとは違うバージョンのSageでオブジェクトを読み込むこともできる場合が多い．
ただし，これは読み書きしたいオブジェクトに関わるコードがバージョン間で大きくは異ならないことが条件となる．
オブジェクトの保存に際しては，その属性の全てがオブジェクトを定義している(ソースコードではなく)クラスと共に保存される．
そのクラスが新バージョンのSageに存在しない場合，配下のオブジェクトを新バージョンでは読み込むことはできない．
しかし古いバージョンで読み込むことはできるはずだから，(``x.__dict__`` で)オブジェクト ``x`` のディクショナリを生成して保存しておけば，それを新しいバージョンで読み込むことができることもある．


テキスト形式で保存する
--------------------------

オブジェクトをASCIIテキスト形式で保存しておくこともできる．
手順は，ファイルを書込みモードで開いて，そこに保存すべきオブジェクトの文字列表現を書き込むだけのことだ(このやり方で複数個のオブジェクトを保存することができる)．
オブジェクトの書込みを終えたら，ファイルをクローズすればよい．


.. skip

::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: f = (x+y)^7
    sage: o = open('file.txt','w')
    sage: o.write(str(f))
    sage: o.close()


.. _section-save:

セッション全体の保存と読み込み
====================================

Sageは，セッション全体を保存し再ロードするための非常に柔軟な機能を備えている．


コマンド ``save_session(セッション名)`` は，現セッション中に定義された全ての変数を、コマンドで指定した ``セッション名`` にディクショナリとして保存する．
(保存を想定していない変数がある場合もまれに見られるが，そうした時はディクショナリに保存されずに終るだけだ．)
保存先は ``.sobj`` ファイルとなり，他の保存済みオブジェクトと全く同じように読み込むことができる．
セッション中に保存したオブジェクトを再びロードすると，変数名をキー，オブジェクトを値とするディクショナリが生成されることになる．


実行中のセッションに ``セッション名`` に定義された変数をロードするには， ``load_session(セッション名)`` コマンドを使う．
このコマンドは現セッションとロードされる側のセッション内容を合併するのであって，現セッションで定義した変数が消去されるわけではない．


まずSageを起動し，変数をいくつか定義しておく．

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: M = ModularSymbols(37)
    sage: a = 389
    sage: t = M.T(2003).matrix(); t.charpoly().factor()
     _4 = (x - 2004) * (x - 12)^2 * (x + 54)^2

次にこのセッションをファイルに保存し，先に定義した変数を残しておく．
``.sobj`` ファイルを確認すると，その大きさは3Kバイトほどとなっている．


.. skip

::

    sage: save_session('misc')
    Saving a
    Saving M
    Saving t
    Saving E
    sage: quit
    was@form:~/tmp$ ls -l misc.sobj
    -rw-r--r--  1 was was 2979 2006-01-28 19:47 misc.sobj

仕上げにSageを再起動し，変数をいくつか追加定義してから，先に保存したセッションを読み込んでみよう．


.. skip

::

    sage: b = 19
    sage: load_session('misc')
    Loading a
    Loading M
    Loading E
    Loading t

保存しておいた変数が再び利用可能になる一方，上で追加した変数 ``b`` は上書きされていないことが分る．


.. skip

::

    sage: M
    Full Modular Symbols space for Gamma_0(37) of weight 2 with sign 0
    and dimension 5 over Rational Field
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational
    Field
    sage: b
    19
    sage: a
    389

