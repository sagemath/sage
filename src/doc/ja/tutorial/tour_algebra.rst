
代数と微積分の基本
==========================

Sageでは，初等的な代数と微積分に関連した多様な演算を実行することができる．
例として，方程式の解を求める，微分や積分を計算する，ラプラス変換の実行などがあげられる．
`Sage Constructions <http://doc.sagemath.org/html/en/constructions/>`_ には，さらに多様な具体例が盛られている．


方程式を解く
-----------------


方程式を解析的に解く
~~~~~~~~~~~~~~~~~~~~~~~~~

``solve``  関数を使って方程式の解を求めることができる．
これを使うには、まず変数を定義し、ついで対象とする方程式(または方程式系)と解くべき変数を ``solve`` の引数として指定する:


::

    sage: x = var('x')
    sage: solve(x^2 + 3*x + 2, x)
    [x == -2, x == -1]

解くべき変数を変更して，解を他の変数で表わすこともできる:


::

    sage: x, b, c = var('x b c')
    sage: solve([x^2 + b*x + c == 0],x)
    [x == -1/2*b - 1/2*sqrt(b^2 - 4*c), x == -1/2*b + 1/2*sqrt(b^2 - 4*c)]


多変数の方程式を解くことも可能だ:

::

    sage: x, y = var('x, y')
    sage: solve([x+y==6, x-y==4], x, y)
    [[x == 5, y == 1]]


次のJason Groutによる例題では、Sageを使って連立非線形方程式を解く．まず，この連立方程式の解を記号的に求めてみよう:


::

    sage: var('x y p q')
    (x, y, p, q)
    sage: eq1 = p+q==9
    sage: eq2 = q*y+p*x==-6
    sage: eq3 = q*y^2+p*x^2==24
    sage: solve([eq1,eq2,eq3,p==1],p,q,x,y)
    [[p == 1, q == 8, x == -4/3*sqrt(10) - 2/3, y == 1/6*sqrt(10) - 2/3], [p == 1, q == 8, x == 4/3*sqrt(10) - 2/3, y == -1/6*sqrt(10) - 2/3]]


解の数値近似を求めるには，やり方を変えて:

.. link

::

    sage: solns = solve([eq1,eq2,eq3,p==1],p,q,x,y, solution_dict=True)
    sage: [[s[p].n(30), s[q].n(30), s[x].n(30), s[y].n(30)] for s in solns]
    [[1.0000000, 8.0000000, -4.8830369, -0.13962039],
     [1.0000000, 8.0000000, 3.5497035, -1.1937129]]


``n`` 関数は解の数値的近似値を表示する. ``n`` の引数は数値精度を表わすビット数を指定している．


方程式を数値的に解く
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

目的の方程式(または方程式系)に対し ``solve`` では厳密解を求めることができないというのは珍しいことではない．
そうした場合には ``find_root`` を使って数値解を求めることができる．
例えば，以下に示す方程式については ``solve`` は何も役に立つことを教えてくれない．

::

    sage: theta = var('theta')
    sage: solve(cos(theta)==sin(theta), theta)
    [sin(theta) == cos(theta)]


しかし代りに ``find_root`` を使えば， :math:`0 < \phi < \pi/2` の範囲で上の方程式の数値解を求めることができる.


::

    sage: phi = var('phi')
    sage: find_root(cos(phi)==sin(phi),0,pi/2)
    0.785398163397448...


微分，積分，その他
----------------------------------

Sageで多様な関数の微分と積分を計算することができる．
例えば :math:`\sin(u)` を :math:`u` で微分するには，以下のようにする:

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

:math:`\sin(x^2)` の4次微分を計算するには:


::

    sage: diff(sin(x^2), x, 4)
    16*x^4*sin(x^2) - 48*x^2*cos(x^2) - 12*sin(x^2)


:math:`x^2+17y^2` の `x` と `y` それぞれによる偏微分を計算するには:


::

    sage: x, y = var('x,y')
    sage: f = x^2 + 17*y^2
    sage: f.diff(x)
    2*x
    sage: f.diff(y)
    34*y


次は不定積分と定積分だ． :math:`\int x\sin(x^2)\, dx` と :math:`\int_0^1 \frac{x}{x^2+1}\, dx` を計算してみよう．


::

    sage: integral(x*sin(x^2), x)
    -1/2*cos(x^2)
    sage: integral(x/(x^2+1), x, 0, 1)
    1/2*log(2)

:math:`\frac{1}{x^2-1}` の部分分数展開を求めるには:

::

    sage: f = 1/((1+x)*(x-1))
    sage: f.partial_fraction(x)
    -1/2/(x + 1) + 1/2/(x - 1)


.. _section-systems:

微分方程式を解く
------------------------------

Sageを使って常微分方程式を研究することもできる． :math:`x'+x-1=0` を解くには:
::

    sage: t = var('t')            # 変数 t を定義
    sage: x = function('x')(t)     # x を t の関数とする
    sage: DE = diff(x, t) + x - 1
    sage: desolve(DE, [x,t])
    (_C + e^t)*e^(-t)


ここでSageはMaxima [Max]_ とインターフェイスしているので，その出力もこれまで見てきたSageの出力とは若干違っている．
上の結果は，上の微分方程式の一般解が :math:`x(t) = e^{-t}(e^{t}+c)` であることを示している．

ラプラス変換を実行することができる．
:math:`t^2e^t -\sin(t)` のラプラス変換は以下のような手順を踏む:

::

    sage: s = var("s")
    sage: t = var("t")
    sage: f = t^2*exp(t) - sin(t)
    sage: f.laplace(t,s)
    -1/(s^2 + 1) + 2/(s - 1)^3


もう少し手間のかかる問題を考えてみよう．
左端が壁に固定された連成バネ各々の、平衡位置からの変位


::

    |------\/\/\/\/\---|mass1|----\/\/\/\/\/----|mass2|
             spring1               spring2

は、連立2階微分方程式


.. math::

    m_1 x_1'' + (k_1+k_2) x_1 - k_2 x_2 = 0

    m_2 x_2''+ k_2 (x_2-x_1) = 0,

でモデル化される．
ここで :math:`m_{i}` はおもり *i* の質量， :math:`x_{i}` はそのおもり *i* の平衡位置からの変位，そして :math:`k_{i}` はバネ *i* のバネ定数である．


**例題:** 上の問題で各パラメータの値を :math:`m_{1}=2`, :math:`m_{2}=1`, :math:`k_{1}=4`, :math:`k_{2}=2`, :math:`x_{1}(0)=3`, :math:`x_{1}'(0)=0`, :math:`x_{2}(0)=3`, :math:`x_{2}'(0)=0` と置き，Sageを使って解いてみよう．

**解法:** まず1番目の方程式をラプラス変換する(記号は :math:`x=x_{1}`, :math:`y=x_{2}` に変える):

::

    sage: de1 = maxima("2*diff(x(t),t, 2) + 6*x(t) - 2*y(t)")
    sage: lde1 = de1.laplace("t","s"); lde1.sage()
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

この出力は読みにくいけれども，意味しているのは

.. math:: -2x'(0) + 2s^2 \cdot X(s) - 2sx(0) - 2Y(s) + 6X(s) = 0

ということだ(ここでは小文字名の関数 :math:`x(t)` のラプラス変換が大文字名の関数 :math:`X(s)` となっている)．
2番目の方程式もラプラス変換してやると:


::

    sage: t,s = SR.var('t,s')
    sage: x = function('x')
    sage: y = function('y')
    sage: f = 2*x(t).diff(t,2) + 6*x(t) - 2*y(t)
    sage: f.laplace(t,s)
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

意味するところは

.. math:: -Y'(0) + s^2Y(s) + 2Y(s) - 2X(s) - sy(0) = 0.

初期条件 :math:`x(0)`, :math:`x'(0)`, :math:`y(0)` ，および :math:`y'(0)` を代入して得られる2つの方程式を `X` と `Y` について解く:

::

    sage: var('s X Y')
    (s, X, Y)
    sage: eqns = [(2*s^2+6)*X-2*Y == 6*s, -2*X +(s^2+2)*Y == 3*s]
    sage: solve(eqns, X,Y)
    [[X == 3*(s^3 + 3*s)/(s^4 + 5*s^2 + 4),
      Y == 3*(s^3 + 5*s)/(s^4 + 5*s^2 + 4)]]

この解の逆ラプラス変換を行なうと:


::

    sage: var('s t')
    (s, t)
    sage: inverse_laplace((3*s^3 + 9*s)/(s^4 + 5*s^2 + 4),s,t)
    cos(2*t) + 2*cos(t)
    sage: inverse_laplace((3*s^3 + 15*s)/(s^4 + 5*s^2 + 4),s,t)
    -cos(2*t) + 4*cos(t)


というわけで，求めていた解は

.. math:: x_1(t) = \cos(2t) + 2\cos(t), \quad x_2(t) = 4\cos(t) - \cos(2t).

これを媒介変数プロットするには

::

    sage: t = var('t')
    sage: P = parametric_plot((cos(2*t) + 2*cos(t), 4*cos(t) - cos(2*t) ),
    ....: (t, 0, 2*pi), rgbcolor=hue(0.9))
    sage: show(P)

各成分ごとにプロットするには


::

    sage: t = var('t')
    sage: p1 = plot(cos(2*t) + 2*cos(t), (t,0, 2*pi), rgbcolor=hue(0.3))
    sage: p2 = plot(4*cos(t) - cos(2*t), (t,0, 2*pi), rgbcolor=hue(0.6))
    sage: show(p1 + p2)


プロットについては :ref:`section-plot` 節の，もう少し詳しい説明を見てほしい．
微分方程式については [NagleEtAl2004]_ の5.5節にもっと詳しい解説がある．


オイラーによる連立微分方程式の解法
----------------------------------------------------

次の例では，1階および2階微分方程式に対するオイラーの解法を具体的に解説する．
手始めに1階微分方程式に対する解法の基本的アイデアを復習しておこう．初期値問題が

.. math::

    y'=f(x,y), \quad y(a)=c,

のような形式で与えられており， :math:`b>a` を満足する :math:`x=b` における解の近似値を求めたいものとする．

微分係数の定義から

.. math::  y'(x) \approx \frac{y(x+h)-y(x)}{h},

ここで :math:`h>0` は与えるべき小さな量である．
この近似式と先の微分方程式を組み合わせると :math:`f(x,y(x))\approx \frac{y(x+h)-y(x)}{h}` が得られる．
これを :math:`y(x+h)` について解くと:

.. math::   y(x+h) \approx y(x) + h\cdot f(x,y(x)).


(他にうまい呼び方も思いつかないので) :math:`h \cdot f(x,y(x))` を "補正項" と呼び， :math:`y(x)` を `y` の "更新前項(old)",  :math:`y(x+h)` を `y` の "更新後項(new)"と呼ぶことにすると，上の近似式を

.. math::   y_{new} \approx y_{old} + h\cdot f(x,y_{old}).

と表わすことができる．


ここで `a` から `b` までの区間を `n` ステップに分割すると :math:`h=\frac{b-a}{n}` と書けるから，ここまでの作業から得られた情報を整理して以下の表のようにまとめることができる．


============== =======================   =====================
:math:`x`      :math:`y`                 :math:`h\cdot f(x,y)`
============== =======================   =====================
:math:`a`      :math:`c`                 :math:`h\cdot f(a,c)`
:math:`a+h`    :math:`c+h\cdot f(a,c)`         ...
:math:`a+2h`   ...
...
:math:`b=a+nh` ???                             ...
============== =======================   =====================


我々の目標は，この表の空欄を上から一行づつ全て埋めていき，最終的に :math:`y(b)` のオイラー法による近似である???に到達することである．

連立微分方程式に対する解法もアイデアは似ている．

**例題:** :math:`z''+tz'+z=0`, :math:`z(0)=1`, :math:`z'(0)=0` を満足する :math:`t=1` における :math:`z(t)` を，4ステップのオイラー法を使って数値的に近似してみよう．

ここでは問題の2階常微分方程式を( :math:`x=z`, :math:`y=z'` として)二つの1階微分方程式に分解してからオイラー法を適用することになる。

::

    sage: t,x,y = PolynomialRing(RealField(10),3,"txy").gens()
    sage: f = y; g = -x - y * t
    sage: eulers_method_2x2(f,g, 0, 1, 0, 1/4, 1)
          t                x            h*f(t,x,y)                y       h*g(t,x,y)
          0                1                  0.00                0           -0.25
        1/4              1.0                -0.062            -0.25           -0.23
        1/2             0.94                 -0.12            -0.48           -0.17
        3/4             0.82                 -0.16            -0.66          -0.081
          1             0.65                 -0.18            -0.74           0.022

したがって， :math:`z(1)\approx 0.65` が判る．

点 :math:`(x,y)` をプロットすれば、その曲線としての概形を見ることができる．
それには関数 ``eulers_method_2x2_plot`` を使うが，その前に三つの成分(`t`, `x`, `y`)からなる引数を持つ関数 `f` と `g` を定義しておかなければならない．

::

    sage: f = lambda z: z[2]        # f(t,x,y) = y
    sage: g = lambda z: -sin(z[1])  # g(t,x,y) = -sin(x)
    sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)

この時点で， ``P`` は2系列のプロットを保持していることになる．
`x` と `t` のプロットである ``P[0]`` ， および  `y` と `t` のプロットである ``P[1]`` である．
これら二つをプロットするには、次のようにする:

.. link

::

    sage: show(P[0] + P[1])

(プロットの詳細については :ref:`section-plot` 節を参照．)


特殊関数
-----------------

数種類の直交多項式と特殊関数が，PARI [GAP]_ およびMaxima [Max]_ を援用して実装されている．
詳細についてはSageレファレンスマニュアルの“Orthogonal polynomials"(直交多項式)と“Special functions"(特殊関数)を参照してほしい．

::

    sage: x = polygen(QQ, 'x')
    sage: chebyshev_U(2,x)
    4*x^2 - 1
    sage: bessel_I(1,1).n(250)
    0.56515910399248502720769602760986330732889962162109200948029448947925564096
    sage: bessel_I(1,1).n()
    0.565159103992485
    sage: bessel_I(2,1.1).n()
    0.167089499251049


ここで注意したいのは，Sageではこれらの関数群が専ら数値計算に便利なようにラップ(wrap)されている点だ．
記号処理をする場合には，以下の例のようにMaximaインターフェイスをじかに呼び出してほしい．

::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'
