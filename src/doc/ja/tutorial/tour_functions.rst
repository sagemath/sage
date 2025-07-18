.. _section-functions-issues:

関数まわりの注意点
=================================

関数の定義については紛らわしい側面があって，微積分やプロットなどを行なう際に問題になることがある．
この節で，関連する諸問題について検討してみたい．

Sageで「関数」と呼ばれるべきものを定義する方法は何通りもある:

1. :ref:`section-functions` 節で解説されている方法で，Python関数を定義する．
こうして定義された関数はプロット可能だが，微分積分演算はできない．

::

       sage: def f(z): return z^2
       sage: type(f)
       <... 'function'>
       sage: f(3)
       9
       sage: plot(f, 0, 2)
       Graphics object consisting of 1 graphics primitive


最終行の書法に注目していただきたい．
これを ``plot(f(z), 0, 2)`` としていたら，エラーになっていたはずである．
``z`` は ``f`` 定義におけるダミー変数であって，定義ブロックの外では未定義になるからだ．
むろん ``f(z)`` のみを実行してもエラーになる．
以下のようにすると切り抜けられるが，どんな場合でも通用するとは限らないので要注意だ(下の第4項を参照)．

.. link

::

       sage: var('z')   # zを変数として定義
       z
       sage: f(z)
       z^2
       sage: plot(f(z), 0, 2)
       Graphics object consisting of 1 graphics primitive

こうすると ``f(z)`` はシンボリック表現になる．シンボリック表現については，次の項目で解説する．


2. 「呼び出し可能シンボリック表現」(callable symbolic expression)を定義する．
これはプロットおよび微分積分演算が可能である．

::

       sage: g(x) = x^2
       sage: g        # gはxをx^2に送る
       x |--> x^2
       sage: g(3)
       9
       sage: Dg = g.derivative(); Dg
       x |--> 2*x
       sage: Dg(3)
       6
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(g, 0, 2)
       Graphics object consisting of 1 graphics primitive

``g`` は呼び出し可能シンボリック表現だが， ``g(x)`` の方はこれに関係はあっても異なる種類のオブジェクトである．
やはりプロットと微積分などが可能なのだが，違っている点もあるので注意を要する．
以下の第5項で具体的に説明する．

.. link

::

       sage: g(x)
       x^2
       sage: type(g(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: g(x).derivative()
       2*x
       sage: plot(g(x), 0, 2)
       Graphics object consisting of 1 graphics primitive


3. Sageで定義済みの「初等関数」(calculus function)を使う.
これらはプロット可能で，ちょっと工夫すると微分積分もできるようになる．


::

       sage: type(sin)
       <class 'sage.functions.trig.Function_sin'>
       sage: plot(sin, 0, 2)
       Graphics object consisting of 1 graphics primitive
       sage: type(sin(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(sin(x), 0, 2)
       Graphics object consisting of 1 graphics primitive

そのままでは ``sin`` は微分演算を受けつけない．
少なくとも ``cos`` にはならない．


::

       sage: f = sin
       sage: f.derivative()
       Traceback (most recent call last):
       ...
       AttributeError: ...


``sin`` そのままではなく ``f = sin(x)`` とすると微積分を受けつけるようになるが， もっと手堅いのは ``f(x) = sin(x)`` として呼び出し可能シンボリック表現を定義することである．


::

       sage: S(x) = sin(x)
       sage: S.derivative()
       x |--> cos(x)


まだ注意を要する点が残っているので，説明しておこう:

4. 意図しない評価が起きることがある．

::

       sage: def h(x):
       ....:     if x < 2:
       ....:         return 0
       ....:     else:
       ....:         return x - 2


ここで ``plot(h(x), 0, 4)`` を実行すると，プロットされるのは `y=x-2` で，複数行にわたって定義しておいた ``h`` ではない．
原因を考えてみよう．
コマンド ``plot(h(x), 0, 4)`` が実行されると，まず ``h(x)`` が評価されるが， これは ``x`` が関数 ``h(x)`` に突っ込まれ ``x<2`` が評価されることを意味する．

.. link

::

       sage: type(x<2)
       <class 'sage.symbolic.expression.Expression'>


シンボリック式が評価される際， ``h`` の定義の場合と同じように，その式が明らかに真でないかぎり戻り値は偽になる．
したがって ``h(x)`` は ``x-2`` と評価され，プロットされるのも ``x-2`` になるわけである．


解決策はというと， ``plot(h(x), 0, 4)`` ではなく


.. link


::

       sage: plot(h, 0, 4)
       Graphics object consisting of 1 graphics primitive

を実行せよ，ということになる．


5. 意図せず関数が定数になってしまう．
::

       sage: f = x
       sage: g = f.derivative()
       sage: g
       1


問題は，例えば ``g(3)`` などと実行するとエラーになって， "ValueError: the number of arguments must be less than or equal to 0."と文句をつけてくることだ．

.. link

::

       sage: type(f)
       <class 'sage.symbolic.expression.Expression'>
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>


``g`` は関数ではなく定数になっているので，変数を持たないから何も値を受けつけない．


解決策は何通りかある．

- ``f`` を最初にシンボリック表式として定義しておく．

::

         sage: f(x) = x        #  'f = x'とはしない
         sage: g = f.derivative()
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>


- または ``f`` の定義は元のまま ``g`` をシンボリック表式として定義する．

::

         sage: f = x
         sage: g(x) = f.derivative()  # 'g = f.derivative()'とするかわり
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>


- または ``f`` と ``g`` の定義は元のまま，代入すべき変数を特定する．

::

         sage: f = x
         sage: g = f.derivative()
         sage: g
         1
         sage: g(x=3)    # たんに'g(3)'とはしない
         1


おしまいになったが， ``f = x`` と ``f(x) = x`` 各々に対する微分の相違点を示す方法がまだあった．


::

       sage: f(x) = x
       sage: g = f.derivative()
       sage: g.variables()  # gに属する変数は?
       ()
       sage: g.arguments()  # gに値を送り込むための引数は?
       (x,)
       sage: f = x
       sage: h = f.derivative()
       sage: h.variables()
       ()
       sage: h.arguments()
       ()


ここの例から判るように， ``h(3)`` がエラーになるのは，そもそも ``h`` が引数を受けつけないためである．
