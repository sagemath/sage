.. sage-doctest: needs sage.plot sage.symbolic

.. _section-plot:

プロットする
==================

Sageを使って2次元および3次元のグラフを作成することができる．


2次元プロット
---------------------

Sageの2次元プロット機能を使うと，円，直線，多辺形の描画はもちろん，
直交座標系における関数のプロット，極座標プロット、等高線プロット，ベクトル場プロットを行うことができる．
以下では、その具体例を見ていくことにしよう．さらにSageによるプロットの具体例を見たければ， :ref:`section-systems` 節と :ref:`section-maxima` 節，および `Sage Constructions <http://doc.sagemath.org/html/en/constructions/>`_ を参照してほしい．

次のコマンドは原点を中心とした半径1の黄色い円を描く:
::

    sage: circle((0,0), 1, rgbcolor=(1,1,0))
    Graphics object consisting of 1 graphics primitive

円を塗りつぶすこともできる：

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0), fill=True)
    Graphics object consisting of 1 graphics primitive

円を生成して，それを変数に収めておくこともできる．
ただし，それだけでは描画は実行されない．

::

    sage: c = circle((0,0), 1, rgbcolor=(1,1,0))

描画するには，以下のように ``c.show()`` または ``show(c)`` などとする:

.. link

::

    sage: c.show()

かわりに ``c.save('filename.png')`` などとすれば，プロット ``c`` は画像としてファイルに保存される.

ところで「円」がどちらかと言えば「楕円」に見えるのは，軸が縦横で異なってスケールされるからだ．
これを直すには:

.. link

::

    sage: c.show(aspect_ratio=1)


同じことはコマンド ``show(c, aspect_ratio=1)`` としても可能だし，画像ファイルとして保存したければ ``c.save('filename.png', aspect_ratio=1)`` と実行してもよい．


初等関数のプロットも簡単だ:

::

    sage: plot(cos, (-5,5))
    Graphics object consisting of 1 graphics primitive

変数を特定しておけば，媒介変数プロットも可能になる:

::

    sage: x = var('x')
    sage: parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    Graphics object consisting of 1 graphics primitive

プロットにおける座標軸の交点は，それがグラフの描画範囲にない限り表示されないことに注意しておいてほしい．
描画範囲として十分大きな値を指定するために，科学記法(指数記法)を用いる必要があるかもしれない．

::

    sage: plot(x^2,(x,300,500))
    Graphics object consisting of 1 graphics primitive

複数のプロットをまとめて描画するには，それらを足し合わせればよい:

::

    sage: x = var('x')
    sage: p1 = parametric_plot((cos(x),sin(x)),(x,0,2*pi),rgbcolor=hue(0.2))
    sage: p2 = parametric_plot((cos(x),sin(x)^2),(x,0,2*pi),rgbcolor=hue(0.4))
    sage: p3 = parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    sage: show(p1+p2+p3, axes=false)


塗りつぶし図形を生成するには，まず図形の頂点からなるリストを生成し(以下の例の ``L``)，次に ``polygon`` コマンドで頂点をつないで閉領域を作るとよい．
ここでは，例として緑色の三角形を描画してみる:

::

    sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),
    ....: 2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,3/4,1/2))
    sage: p
    Graphics object consisting of 1 graphics primitive

``show(p, axes=false)`` と入力して座標軸なしの図を表示してみよう．

図形プロットに文字列を加えることができる:

::

    sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),
    ....: 6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))
    sage: t = text("hypotrochoid", (5,4), rgbcolor=(1,0,0))
    sage: show(p+t)


次に出てくるarcsinのグラフは、微積分の教師が黒板にたびたび描くものだ．
主値だけではなく複数の分岐を含めてプロットするには， :math:`x` が :math:`-2\pi` から :math:`2\pi` までの :math:`y=\sin(x)` のグラフを傾き45度の線について反転させて描いてやればよい．
これをSageで実行するには，以下のコマンドを使う:

::

    sage: v = [(sin(x),x) for x in srange(-2*float(pi),2*float(pi),0.1)]
    sage: line(v)
    Graphics object consisting of 1 graphics primitive

正接(tan)関数はsin関数よりも値域が広いので，arcsinと同じ作戦で逆正接関数のグラフを描くには *x* -軸の最大値と最小値を調節してやる必要がある:

::

    sage: v = [(tan(x),x) for x in srange(-2*float(pi),2*float(pi),0.01)]
    sage: show(line(v), xmin=-20, xmax=20)

Sageでは、(対象となる関数は限られるが)極座標プロット、等高線プロット、ベクトル場プロットも可能だ．
ここでは，例として等高線プロットを見ておこう:

::

    sage: f = lambda x,y: cos(x*y)
    sage: contour_plot(f, (-4, 4), (-4, 4))
    Graphics object consisting of 1 graphics primitive


3次元プロット
-----------------------

Sageでは3次元プロットも作成することができる．
ノートブック上でもREPL(コマンドライン)上でも，3次元プロットの表示はデフォルトでオープンソースパッケージ [ThreeJS]_ によって行なわれる．
Jmolではマウスによる描画の回転と拡大縮小が可能だ．

``plot3d`` を使って `f(x, y) = z` 形式の関数をプロットしてみよう:

::

    sage: x, y = var('x,y')
    sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
    Graphics3d Object

代りに ``parametric_plot3d`` を使い， `x, y, z` 各々が1あるいは2個のパラメター(いわゆる媒介変数，記号 `u` や `v` などが使われることが多い)で決定されるパラメトリック曲面として描画することもできる．
上の関数を媒介変数表示してプロットするには:

::

    sage: u, v = var('u, v')
    sage: f_x(u, v) = u
    sage: f_y(u, v) = v
    sage: f_z(u, v) = u^2 + v^2
    sage: parametric_plot3d([f_x, f_y, f_z], (u, -2, 2), (v, -2, 2))
    Graphics3d Object

Sageで3次元曲面プロットを行うための第三の方法が ``implicit_plot3d`` の使用で，これは(空間内の点の集合を定義する) `f(x, y, z) = 0` を満足する関数の等高線を描画する．
ここでは古典的な表式を使って球面を作画してみよう:

::

    sage: x, y, z = var('x, y, z')
    sage: implicit_plot3d(x^2 + y^2 + z^2 - 4, (x,-2, 2), (y,-2, 2), (z,-2, 2))
    Graphics3d Object

以下で，さらにいくつかの3次元プロットを示しておこう:

.. `Yellow Whitney's umbrella <http://en.wikipedia.org/wiki/Whitney_umbrella>`__:

`ホィットニーの傘 <http://en.wikipedia.org/wiki/Whitney_umbrella>`__:

::

    sage: u, v = var('u,v')
    sage: fx = u*v
    sage: fy = u
    sage: fz = v^2
    sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -1, 1),
    ....: frame=False, color="yellow")
    Graphics3d Object

`クロスキャップ(十字帽) <http://en.wikipedia.org/wiki/Cross-cap>`__:

::

    sage: u, v = var('u,v')
    sage: fx = (1+cos(v))*cos(u)
    sage: fy = (1+cos(v))*sin(u)
    sage: fz = -tanh((2/3)*(u-pi))*sin(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....: frame=False, color="red")
    Graphics3d Object

ねじれトーラス(twisted torus):

::

    sage: u, v = var('u,v')
    sage: fx = (3+sin(v)+cos(u))*cos(2*v)
    sage: fy = (3+sin(v)+cos(u))*sin(2*v)
    sage: fz = sin(u)+2*cos(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....: frame=False, color="red")
    Graphics3d Object

レムニスケート(連珠形, lemniscate):

::

    sage: x, y, z = var('x,y,z')
    sage: f(x, y, z) = 4*x^2 * (x^2 + y^2 + z^2 + z) + y^2 * (y^2 + z^2 - 1)
    sage: implicit_plot3d(f, (x, -0.5, 0.5), (y, -1, 1), (z, -1, 1))
    Graphics3d Object
