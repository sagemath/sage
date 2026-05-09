********
代数几何
********

.. index::
   pair: elliptic curve; point counting

曲线上点的计数
==============

如何在 Sage 中计算有限域上椭圆曲线上点的数量？

对于素有限域，有小步大步法 (baby step giant step method) 和 SEA (Schoof-Elkies-Atkin) 算法
（由 Christophe Doche 和 Sylvain Duquesne 在 PARI 中实现）。
以下是从参考手册中摘取的示例：

::

    sage: E = EllipticCurve(GF(10007),[1,2,3,4,5])
    sage: E.cardinality()
    10076

``E.points()`` 命令将返回有理点的实际列表。

如何在有限域上计算平面曲线上有理点的数量？
``rational_points`` 命令使用简单枚举算法生成有理点。以下是语法示例：

::

    sage: x,y,z = PolynomialRing(GF(5), 3, 'xyz').gens()
    sage: C = Curve(y^2*z^7 - x^9 - x*z^8); C
    Projective Plane Curve over Finite Field of size 5 defined by -x^9 + y^2*z^7 - x*z^8
    sage: C.rational_points()
    [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1), (3 : 1 : 1), (3 : 4 : 1)]
    sage: C.rational_points(algorithm="bn")
    [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1), (3 : 1 : 1), (3 : 4 : 1)]

选项 ``algorithm="bn"`` 使用 Sage 的 Singular 接口并调用 ``brnoeth`` 包。

下面是将 Sage 的 ``rational_points`` 应用于 `GF(8)` 上的 Klein 四次方程的另一个例子。

::

    sage: x, y, z = PolynomialRing(GF(8,'a'), 3, 'xyz').gens()
    sage: f = x^3*y+y^3*z+x*z^3
    sage: C = Curve(f); C
    Projective Plane Curve over Finite Field in a of size 2^3 defined by x^3*y + y^3*z + x*z^3
    sage: C.rational_points()
    [(0 : 0 : 1),
     (0 : 1 : 0),
     (1 : 0 : 0),
     (1 : a : 1),
     (1 : a^2 : 1),
     (1 : a^2 + a : 1),
     (a : 1 : 1),
     (a : a^2 : 1),
     (a : a^2 + 1 : 1),
     (a + 1 : a + 1 : 1),
     (a + 1 : a^2 : 1),
     (a + 1 : a^2 + a + 1 : 1),
     (a^2 : 1 : 1),
     (a^2 : a^2 + a : 1),
     (a^2 : a^2 + a + 1 : 1),
     (a^2 + 1 : a + 1 : 1),
     (a^2 + 1 : a^2 + 1 : 1),
     (a^2 + 1 : a^2 + a : 1),
     (a^2 + a : 1 : 1),
     (a^2 + a : a : 1),
     (a^2 + a : a + 1 : 1),
     (a^2 + a + 1 : a : 1),
     (a^2 + a + 1 : a^2 + 1 : 1),
     (a^2 + a + 1 : a^2 + a + 1 : 1)]

其他方法
--------


-  对于平面曲线，你可以使用 Singular 的 ``closed_points`` 命令。
   输入是 `2` 变量环 `F[x,y]` 中曲线 `X` 的消失理想 `I`。
   ``closed_points`` 命令返回一个素理想列表（每个都是 Gröbner 基），
   对应于 `V(I)` 的（不同仿射闭合）点。以下是示例：

   .. skip

   ::

       sage: singular_console()
                            SINGULAR                             /  Development
        A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                              0<
            by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
       FB Mathematik der Universitaet, D-67653 Kaiserslautern    \
       // ** executing /home/wdj/sagefiles/sage-0.9.4/local/LIB/.singularrc
       > LIB "brnoeth.lib";
       > ring s = 2,(x,y),lp;
       > ideal I = x4+x,y4+y;
       > list L = closed_points(I);
       > L;
       [1]:
          _[1] = y
          _[2] = x
       [2]:
          _[1] = y
          _[2] = x+1
       [3]:
          _[1] = y
          _[2] = x2+x+1
       [4]:
          _[1] = y+1
          _[2] = x
       [5]:
          _[1] = y+1
          _[2] = x+1
       [6]:
          _[1] = y+1
          _[2] = x2+x+1
       [7]:
          _[1] = y2+y+1
          _[2] = x+1
       [8]:
          _[1] = y2+y+1
          _[2] = x
       [9]:
          _[1] = y2+y+1
          _[2] = x+y
       [10]:
          _[1] = y2+y+1
          _[2] = x+y+1
       > Auf Wiedersehen.

   ::

       sage: singular.lib("brnoeth.lib")
       sage: s = singular.ring(2,'(x,y)','lp')
       sage: I = singular.ideal('x^4+x', 'y^4+y')
       sage: L = singular.closed_points(I)
       sage: # Here you have all the points :
       sage: L       # random
       [1]:
          _[1]=y+1
          _[2]=x+1
       ...
       sage: l=[L[k].sage() for k in [1..10]]; len(l) # there are 10 points
       10
       sage: r=sorted(l[0].ring().gens()); r
       [y, x]
       sage: r in [t.gens() for t in l] #  one of them is given by [y,x]
       True

-  另一种计算有理点的方法是使用 Singular 的 ``NSplaces`` 命令。
   以下是用该方法在 `GF(8)` 上计算 Klein 四次方程的示例：

   ::

       sage: singular.LIB("brnoeth.lib")
       sage: s = singular.ring(2,'(x,y)','lp')
       sage: f = singular.poly('x3y+y3+x')
       sage: klein1 = f.Adj_div(); print(klein1)
       [1]:
          [1]:
             // coefficients: ZZ/2...
       // number of vars : 2
       //        block   1 : ordering lp
       //                  : names    x y
       //        block   2 : ordering C
       ...
       sage: # define a curve X = {f = 0} over GF(2)
       sage: klein2 = singular.NSplaces(3,klein1)
       sage: print(singular.eval('extcurve(3,%s)'%klein2.name()))
       Total number of rational places : NrRatPl = 23
       ...
       sage: klein3 = singular.extcurve(3, klein2)

   上面我们在 Singular 中定义了一条 `GF(8)` 上的曲线 `X = \{f = 0\}`。

   .. link

   ::

       sage: print(klein1)
       [1]:
          [1]:
             // coefficients: ZZ/2...
       // number of vars : 2
       //        block   1 : ordering lp
       //                  : names    x y
       //        block   2 : ordering C
          [2]:
             // coefficients: ZZ/2...
       // number of vars : 3
       //        block   1 : ordering lp
       //                  : names    x y z
       //        block   2 : ordering C
       [2]:
          4,3
       [3]:
          [1]:
             1,1
          [2]:
             1,2
       [4]:
          0
       [5]:
          [1]:
             [1]:
                // coefficients: ZZ/2...
       // number of vars : 3
       //        block   1 : ordering ls
       //                  : names    x y t
       //        block   2 : ordering C
             [2]:
                1,1
       sage: print(klein1[3])
       [1]:
          1,1
       [2]:
          1,2

   对于度为 `3` 的地方：

   .. link

   ::

       sage: print(klein2[3])
       [1]:
          1,1
       [2]:
          1,2
       [3]:
          3,1
       [4]:
          3,2
       [5]:
          3,3
       [6]:
          3,4
       [7]:
          3,5
       [8]:
          3,6
       [9]:
          3,7

   下面的每个点都是成对的：（度，点索引号）。

   .. link

   ::

       sage: print(klein3[3])
       [1]:
          1,1
       [2]:
          1,2
       [3]:
          3,1
       [4]:
          3,2
       [5]:
          3,3
       [6]:
          3,4
       [7]:
          3,5
       [8]:
          3,6
       [9]:
          3,7

   实际获取 `X(GF(8))` 的点：

   .. link

   ::

       sage: R = klein3[1][5]
       sage: R.set_ring()
       sage: singular("POINTS;")
       [1]:
          [1]:
             0
          [2]:
             1
          [3]:
             0
       [2]:
          [1]:
             1
          [2]:
             0
          [3]:
             0
       ...

   加上另外 21 个点（已省略）。总共有 `23` 个有理点。

.. index:: Riemann-Roch space

使用 Singular 计算 Riemann-Roch 空间
====================================

为了计算域 `F` 上曲线上的因子 `D` 的 Riemann-Roch 空间基，
可以使用 Sage 封装的 ``riemann_roch_basis`` 方法，它是 Singular 实现的 Brill Noether 算法。
注意，这个封装当前仅在 `F` 是素数且因子 `D` 在有理点上受支持时才有效。
下面是如何使用 ``riemann_roch_basis`` 的示例，以及如何使用 Singular 本身来帮助理解封装的工作方式。

-  使用 ``riemann_roch_basis``:

   ::

       sage: x, y, z = PolynomialRing(GF(5), 3, 'xyz').gens()
       sage: f = x^7 + y^7 + z^7
       sage: X = Curve(f); pts = X.rational_points()
       sage: D = X.divisor([ (3, pts[0]), (-1,pts[1]), (10, pts[5]) ])
       sage: X.riemann_roch_basis(D)
       [(-2*x + y)/(x + y), (-x + z)/(x + y)]

-  使用 Singular 的 ``BrillNoether`` 命令
   （具体内容请参见 Singular 在线文档的 Brill-Noether 章节
   https://www.singular.uni-kl.de/Manual/4-3-0/sing_2254.htm 和论文{CF}）：

   ::

       sage: singular.LIB('brnoeth.lib')
       sage: _ = singular.ring(5,'(x,y)','lp')
       sage: print(singular.eval("list X = Adj_div(-x5+y2+x);"))
       Computing affine singular points ...
       Computing all points at infinity ...
       Computing affine singular places ...
       Computing singular places at infinity ...
       Computing non-singular places at infinity ...
       Adjunction divisor computed successfully
       <BLANKLINE>
       The genus of the curve is 2
       sage: print(singular.eval("X = NSplaces(1,X);"))
       Computing non-singular affine places of degree 1 ...
       sage: print(singular("X[3];"))
       [1]:
          1,1
       [2]:
          1,2
       [3]:
          1,3
       [4]:
          1,4
       [5]:
          1,5
       [6]:
          1,6

   上述列表中，每个整数对中的第一个整数表示点的度数 `d`。
   第二个整数是该点在环 X[5][`d`][1] 的 POINTS 列表中的索引。
   注意，每次运行算法时，这个列表的顺序都不相同，
   例如上面列表中的 `1`, `1` 每次可能指示不同的有理点。
   通过定义一个与 X[3] 长度相同的整数列表 `G`，可以指定一个因子。
   如果 X[3] 的第 `k` 项为 `d`, `i`，则 `G` 的第 `k` 项表示
   该因子在环 X[5][`d`][1] 的 POINTS 列表中第 `i` 个点上的重数。
   接下来，我们定义一个度为 12 的“随机”因子并计算其 Riemann-Roch 空间基：

   .. link

   ::

       sage: singular.eval("intvec G = 4,4,4,0,0,0;")
       ''
       sage: singular.eval("def R = X[1][2];")
       ''
       sage: singular.eval("setring R;")
       ''
       sage: print(singular.eval("list LG = BrillNoether(G,X);"))
       Forms of degree 6 :
       28
       <BLANKLINE>
       Vector basis successfully computed
       <BLANKLINE>


.. index::
   pair: codes; algebraic-geometric

AG 码
-----

Sage 可以通过调用 Singular 的 BrillNoether 算法计算 Riemann-Roch 空间 `L(D)=L_X(D)` 的基，
从而计算 AG 码 `C=C_X(D,E)`。
除了曲线 `X` 和因子 `D`，还必须指定求值因子 `E`。

请注意，自从 ``riemann_roch_basis`` 封装被修复后，本节尚未更新。
请参阅上文中，了解如何正确定义 Singular 的 ``BrillNoether`` 命令的因子。

这里有一个示例，计算相关 AG 码的生成矩阵。这次我们使用 Singular 的 ``AGCode_L`` 命令：

::

    sage: singular.LIB('brnoeth.lib')
    sage: singular.eval("ring s = 2,(x,y),lp;")
    ''
    sage: print(singular.eval("list HC = Adj_div(x3+y2+y);"))
    Computing affine singular points ...
    Computing all points at infinity ...
    Computing affine singular places ...
    Computing singular places at infinity ...
    Computing non-singular places at infinity ...
    Adjunction divisor computed successfully
    <BLANKLINE>
    The genus of the curve is 1
    sage: print(singular.eval("list HC1 = NSplaces(1..2,HC);"))
    Computing non-singular affine places of degree 1 ...
    Computing non-singular affine places of degree 2 ...
    sage: print(singular.eval("HC = extcurve(2,HC1);"))
    Total number of rational places : NrRatPl = 9

我们将以下内容设置为 ``junk`` 以丢弃输出::

    sage: junk = singular.eval("intvec G = 5;")      # the rational divisor G = 5*HC[3][1]
    sage: junk = singular.eval("def R = HC[1][2];")
    sage: singular.eval("setring R;")
    ''

向量 `G` 表示因子“无穷远点的 5 倍”。

.. index:: Riemann-Roch space

接下来，我们计算 Riemann-Roch 空间。

.. link

::

    sage: print(singular.eval("BrillNoether(G,HC);"))
    Forms of degree 3 :
    10
    <BLANKLINE>
    Vector basis successfully computed
    <BLANKLINE>
    [1]:
       _[1]=x
       _[2]=z
    [2]:
       _[1]=y
       _[2]=z
    [3]:
       _[1]=1
       _[2]=1
    [4]:
       _[1]=y2+yz
       _[2]=xz
    [5]:
       _[1]=y3+y2z
       _[2]=x2z

这是 Riemann-Roch 空间的基，其中每个函数对表示商（第一个函数除以第二个函数）。
每一个基元素都会在特定点进行求值以构建码的生成矩阵。接下来我们构建这些点。

.. skip

::

    sage: singular.eval("def R = HC[1][5];")
    '// ** redefining R **'
    sage: singular.eval("setring R;")
    ''
    sage: print(singular.eval("POINTS;"))
    [1]:
       [1]:
          0
       [2]:
          1
       [3]:
          0
    [2]:
       [1]:
          0
       [2]:
          1
       [3]:
          1
    [3]:
       [1]:
          0
       [2]:
          0
       [3]:
          1
    [4]:
       [1]:
          (a+1)
       [2]:
          (a)
       [3]:
          1
    ...

再加上 `5` 个，曲线上总共有 `9` 个有理点。
我们使用这些点的子集（除第一个点外）定义我们的“求值因子” `D`：

.. skip

::

    sage: singular.eval("def ER = HC[1][4];")
    ''
    sage: singular.eval("setring ER;")
    ''
    sage: # D = sum of the rational places no. 2..9 over F_4
    sage: singular.eval("intvec D = 2..9;")
    ''
    sage: # let us construct the corresponding evaluation AG code :
    sage: print(singular.eval("matrix C = AGcode_L(G,D,HC);"))
    Forms of degree 3 :
    10
    <BLANKLINE>
    Vector basis successfully computed
    <BLANKLINE>
    sage: # here is a linear code of type [8,5,> = 3] over F_4
    sage: print(singular.eval("print(C);"))
    0,0,(a+1),(a),  1,  1,    (a),  (a+1),
    1,0,(a),  (a+1),(a),(a+1),(a),  (a+1),
    1,1,1,    1,    1,  1,    1,    1,
    0,0,(a),  (a+1),1,  1,    (a+1),(a),
    0,0,1,    1,    (a),(a+1),(a+1),(a)

这就是我们最终想要的生成矩阵，其中 ``a`` 表示基域 `GF(2)` 上度为 `2` 的域扩张的生成器。

是否可以对其进行“封装”？
