一些更高级的数学
================

代数几何
--------

在 Sage 中可以定义任意代数簇，但有时复杂的功能仅限于在 `\QQ` 或有限域上的环。
例如，我们可以计算两个仿射平面曲线的并集，然后将曲线恢复成该并集的不可约分量。

::

    sage: x, y = AffineSpace(2, QQ, 'xy').gens()
    sage: C2 = Curve(x^2 + y^2 - 1)
    sage: C3 = Curve(x^3 + y^3 - 1)
    sage: D = C2 + C3
    sage: D
    Affine Plane Curve over Rational Field defined by
       x^5 + x^3*y^2 + x^2*y^3 + y^5 - x^3 - y^3 - x^2 - y^2 + 1
    sage: D.irreducible_components()
    [Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x^2 + y^2 - 1,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x^3 + y^3 - 1]

我们还可以通过相交这两条曲线并计算其不可约分量来找到它们的所有交点。

.. link

::

    sage: V = C2.intersection(C3)
    sage: V.irreducible_components()
    [Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       y - 1,
       x,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       y,
       x - 1,
     Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
       x + y + 2,
       2*y^2 + 4*y + 3]

例如，`(1,0)` 和 `(0,1)` 都在两条曲线上（显而易见），
还有一些（二次）点，它们的 `y` 坐标满足 `2y^2 + 4y + 3=0`。

Sage 可以计算三维射影空间中扭曲三次曲线的环理想：

::

    sage: R.<a,b,c,d> = PolynomialRing(QQ, 4)
    sage: I = ideal(b^2-a*c, c^2-b*d, a*d-b*c)
    sage: F = I.groebner_fan(); F
    Groebner fan of the ideal:
    Ideal (b^2 - a*c, c^2 - b*d, -b*c + a*d) of Multivariate Polynomial Ring
    in a, b, c, d over Rational Field
    sage: F.reduced_groebner_bases ()
    [[-c^2 + b*d, -b*c + a*d, -b^2 + a*c],
     [-b*c + a*d, -c^2 + b*d, b^2 - a*c],
     [-c^3 + a*d^2, -c^2 + b*d, b*c - a*d, b^2 - a*c],
     [-c^2 + b*d, b^2 - a*c, b*c - a*d, c^3 - a*d^2],
     [-b*c + a*d, -b^2 + a*c, c^2 - b*d],
     [-b^3 + a^2*d, -b^2 + a*c, c^2 - b*d, b*c - a*d],
     [-b^2 + a*c, c^2 - b*d, b*c - a*d, b^3 - a^2*d],
     [c^2 - b*d, b*c - a*d, b^2 - a*c]]
    sage: F.polyhedralfan()
    Polyhedral fan in 4 dimensions of dimension 4

椭圆曲线
--------

Sage 的椭圆曲线功能包括 PARI 的大部分椭圆曲线功能、
访问 Cremona 在线表中的数据（需要可选数据库包）、
mwrank 功能（即计算全 Mordell-Weil 群的 2 次下降）、
SEA 算法、计算所有同源、许多关于 `\QQ` 上曲线的新代码，
以及 Denis Simon 的一些代数下降软件。

创建椭圆曲线的命令 ``EllipticCurve`` 有多种形式：


-  EllipticCurve([`a_1`, `a_2`, `a_3`, `a_4`, `a_6`]):
   返回如下椭圆曲线

   .. MATH::  y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6,


   其中 `a_i` 被转换为 `a_1` 的父结构。
   如果所有 `a_i` 的父结构都是 `\ZZ`，它们将被转换为 `\QQ`。

-  EllipticCurve([`a_4`, `a_6`]): 与上面相同，但
   `a_1=a_2=a_3=0`。

-  EllipticCurve(label): 返回来自 Cremona 数据库的椭圆曲线，使用给定的（新的！）Cremona 标签。
   标签是一个字符串，例如 ``"11a"`` 或 ``"37b2"``。字母必须是小写（以区分旧标签）。

-  EllipticCurve(j): 返回具有 `j`-不变量 `j` 的椭圆曲线。

-  EllipticCurve(R,
   [`a_1`, `a_2`, `a_3`, `a_4`, `a_6`]):
   创建定义在环 `R` 上的椭圆曲线，给定的 `a_i` 同上。


我们将展示每一个构造函数：

::

    sage: EllipticCurve([0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve([GF(5)(0),0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

    sage: EllipticCurve([1,2])
    Elliptic Curve defined by y^2  = x^3 + x + 2 over Rational Field

    sage: EllipticCurve('37a')
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve_from_j(1)
    Elliptic Curve defined by y^2 + x*y = x^3 + 36*x + 3455 over Rational Field

    sage: EllipticCurve(GF(5), [0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

点 `(0,0)` 是椭圆曲线 `E` 上的一点，定义为 `y^2 + y = x^3 - x`。
要在 Sage 中创建该点，输入 ``E([0,0])``。Sage 可以在此椭圆曲线上添加点
（椭圆曲线支持一个加法群结构，其中无穷远点为零元素，曲线上三个共线点之和为零）：

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
    sage: P = E([0,0])
    sage: P + P
    (1 : 0 : 1)
    sage: 10*P
    (161/16 : -2065/64 : 1)
    sage: 20*P
    (683916417/264517696 : -18784454671297/4302115807744 : 1)
    sage: E.conductor()
    37

复数域上的椭圆曲线由 `j`-不变量参数化。Sage 计算 `j`-不变量如下：

::

    sage: E = EllipticCurve([0,0,0,-4,2]); E
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: E.conductor()
    2368
    sage: E.j_invariant()
    110592/37

如果我们创建一个具有与 `E` 相同 `j`-不变量的曲线，它不一定与 `E` 同构。
在以下示例中，这些曲线不相同，因为它们的导数不同。

::

    sage: F = EllipticCurve_from_j(110592/37)
    sage: F.conductor()
    37

然而，通过对 `F` 进行 2 次扭转可以得到一个与其同构的曲线。

.. link

::

    sage: G = F.quadratic_twist(2); G
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: G.conductor()
    2368
    sage: G.j_invariant()
    110592/37

我们可以计算椭圆曲线的 `L`-级数或模形式 `\sum_{n=0}^\infty a_nq^n` 的系数 `a_n`。
此计算使用 PARI C 库：

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E.anlist(30)
    [0, 1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0, -4,
     3, 10, 2, 0, -1, 4, -9, -2, 6, -12]
    sage: v = E.anlist(10000)

对于 `n\leq 10^5`，计算所有 `a_n` 仅需几秒：

.. skip

::

    sage: %time v = E.anlist(100000)
    CPU times: user 0.98 s, sys: 0.06 s, total: 1.04 s
    Wall time: 1.06

椭圆曲线可以使用它们的 Cremona 标签构造。
这会预加载椭圆曲线的秩、Tamagawa 数、调节器等信息。

::

    sage: E = EllipticCurve("37b2")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational
    Field
    sage: E = EllipticCurve("389a")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x  over Rational Field
    sage: E.rank()
    2
    sage: E = EllipticCurve("5077a")
    sage: E.rank()
    3

我们也可以直接访问 Cremona 数据库。

::

    sage: db = sage.databases.cremona.CremonaDatabase()
    sage: db.curves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1], 'b1': [[0, 1, 1, -23, -50], 0, 3]}
    sage: db.allcurves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1],
     'b1': [[0, 1, 1, -23, -50], 0, 3],
     'b2': [[0, 1, 1, -1873, -31833], 0, 1],
     'b3': [[0, 1, 1, -3, 1], 0, 3]}

从数据库返回的对象不是 ``EllipticCurve`` 类型。
它们是数据库中的元素，只有几个字段而已。
Cremona 数据库有一个小型版本，默认随 Sage 一起分发，包含有关导子(conductor) `\leq 10000` 的椭圆曲线的有限信息。
还有一个大型可选版本，包含有关所有导子不超过 `120000` 的曲线的大量数据（截至 2005 年 10 月）。
Sage 还有一个巨大的（2GB）可选数据库包，包含 Stein-Watkins 数据库中数亿条椭圆曲线数据。

狄利克雷特征
------------

*Dirichlet 特征* 是同态 `(\ZZ/N\ZZ)^* \to R^*` 的扩展，
对于某个环 `R`，可以通过将满足 `\gcd(N,x)>1` 的整数 `x` 映射到 0
从而得到一个 `\ZZ \to R` 的映射。

::

    sage: G = DirichletGroup(12)
    sage: G.list()
    [Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1,
    Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1]
    sage: G.gens()
    (Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1)
    sage: len(G)
    4

创建该群之后，我们继续创建一个元素并进行计算。

.. link

::

    sage: G = DirichletGroup(21)
    sage: chi = G.1; chi
    Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6
    sage: chi.values()
    [0, 1, zeta6 - 1, 0, -zeta6, -zeta6 + 1, 0, 0, 1, 0, zeta6, -zeta6, 0, -1,
     0, 0, zeta6 - 1, zeta6, 0, -zeta6 + 1, -1]
    sage: chi.conductor()
    7
    sage: chi.modulus()
    21
    sage: chi.order()
    6
    sage: chi(19)
    -zeta6 + 1
    sage: chi(40)
    -zeta6 + 1

还可以计算伽罗瓦群 `\text{Gal}(\QQ(\zeta_N)/\QQ)` 对这些特征的作用，
以及对应于模数分解的直积分解。

.. link

::

    sage: chi.galois_orbit()
    [Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> -zeta6 + 1,
     Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6]

    sage: go = G.galois_orbits()
    sage: [len(orbit) for orbit in go]
    [1, 2, 2, 1, 1, 2, 2, 1]

    sage: G.decomposition()
    [Group of Dirichlet characters modulo 3 with values in Cyclotomic Field of order 6 and degree 2,
     Group of Dirichlet characters modulo 7 with values in Cyclotomic Field of order 6 and degree 2]

接下来，我们构造模 20 的狄利克雷特征群，但其值在 `\QQ(i)` 中：

::

    sage: K.<i> = NumberField(x^2+1)
    sage: G = DirichletGroup(20,K)
    sage: G
    Group of Dirichlet characters modulo 20 with values in Number Field in i with defining polynomial x^2 + 1


接下来我们计算 ``G`` 的几个不变量：

.. link

::

    sage: G.gens()
    (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1,
    Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> i)

    sage: G.unit_gens()
    (11, 17)
    sage: G.zeta()
    i
    sage: G.zeta_order()
    4

下面这个例子中，我们创建了一个值在数域中的狄利克雷特征。通过 ``DirichletGroup`` 的第三个参数明确指定了选择的单位根。

::

    sage: x = polygen(QQ, 'x')
    sage: K = NumberField(x^4 + 1, 'a'); a = K.0
    sage: b = K.gen(); a == b
    True
    sage: K
    Number Field in a with defining polynomial x^4 + 1
    sage: G = DirichletGroup(5, K, a); G
    Group of Dirichlet characters modulo 5 with values in the group of order 8 generated by a in Number Field in a with defining polynomial x^4 + 1
    sage: chi = G.0; chi
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> a^2
    sage: [(chi^i)(2) for i in range(4)]
    [1, a^2, -1, -a^2]

这里 ``NumberField(x^4 + 1, 'a')`` 告诉 Sage 在打印 ``K`` 时使用符号 "a"
（一个定义多项式 `x^4 + 1` 的数域）。此时名称 "a" 尚未声明。
一旦执行 ``a = K.0`` （或等价的 ``a = K.gen()``），符号 "a" 就代表生成多项式 `x^4+1` 的一个根。

模形式
------

Sage 可以进行一些与模形式相关的计算，包括计算维度、模符号空间、Hecke 算子和分解。

有几个函数可以用来计算模形式空间的维度。例如，

::

    sage: from sage.modular.dims import dimension_cusp_forms
    sage: dimension_cusp_forms(Gamma0(11),2)
    1
    sage: dimension_cusp_forms(Gamma0(1),12)
    1
    sage: dimension_cusp_forms(Gamma1(389),2)
    6112

接下来我们展示如何在权重 `12` 和级别 `1` 的模符号空间上计算 Hecke 算子。

::

    sage: M = ModularSymbols(1,12)
    sage: M.basis()
    ([X^8*Y^2,(0,0)], [X^9*Y,(0,0)], [X^10,(0,0)])
    sage: t2 = M.T(2)
    sage: t2
    Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1)
    of weight 12 with sign 0 over Rational Field
    sage: t2.matrix()
    [ -24    0    0]
    [   0  -24    0]
    [4860    0 2049]
    sage: f = t2.charpoly('x'); f
    x^3 - 2001*x^2 - 97776*x - 1180224
    sage: factor(f)
    (x - 2049) * (x + 24)^2
    sage: M.T(11).charpoly('x').factor()
    (x - 285311670612) * (x - 534612)^2

我们还可以创建 `\Gamma_0(N)` 和 `\Gamma_1(N)` 的模符号空间。

::

    sage: ModularSymbols(11,2)
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
     0 over Rational Field
    sage: ModularSymbols(Gamma1(11),2)
    Modular Symbols space of dimension 11 for Gamma_1(11) of weight 2 with
    sign 0 over Rational Field

让我们计算一些特征多项式和 `q` 展开式。

::

    sage: M = ModularSymbols(Gamma1(11),2)
    sage: M.T(2).charpoly('x')
    x^11 - 8*x^10 + 20*x^9 + 10*x^8 - 145*x^7 + 229*x^6 + 58*x^5 - 360*x^4
         + 70*x^3 - 515*x^2 + 1804*x - 1452
    sage: M.T(2).charpoly('x').factor()
    (x - 3) * (x + 2)^2 * (x^4 - 7*x^3 + 19*x^2 - 23*x + 11)
            * (x^4 - 2*x^3 + 4*x^2 + 2*x + 11)
    sage: S = M.cuspidal_submodule()
    sage: S.T(2).matrix()
    [-2  0]
    [ 0 -2]
    sage: S.q_expansion_basis(10)
    [q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)]

我们甚至可以计算带有特征的模符号空间。

::

    sage: G = DirichletGroup(13)
    sage: e = G.0^2
    sage: M = ModularSymbols(e,2); M
    Modular Symbols space of dimension 4 and level 13, weight 2, character
    [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
    sage: M.T(2).charpoly('x').factor()
    (x - zeta6 - 2) * (x - 2*zeta6 - 1) * (x + zeta6 + 1)^2
    sage: S = M.cuspidal_submodule(); S
    Modular Symbols subspace of dimension 2 of Modular Symbols space of
    dimension 4 and level 13, weight 2, character [zeta6], sign 0, over
    Cyclotomic Field of order 6 and degree 2
    sage: S.T(2).charpoly('x').factor()
    (x + zeta6 + 1)^2
    sage: S.q_expansion_basis(10)
    [q + (-zeta6 - 1)*q^2 + (2*zeta6 - 2)*q^3 + zeta6*q^4 + (-2*zeta6 + 1)*q^5 + (-2*zeta6 + 4)*q^6 + (2*zeta6 - 1)*q^8 - zeta6*q^9 + O(q^10)]

以下是 Sage 如何计算 Hecke 算子在模形式空间上的作用的另一个例子。

::

    sage: T = ModularForms(Gamma0(11),2)
    sage: T
    Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of
    weight 2 over Rational Field
    sage: T.degree()
    2
    sage: T.level()
    11
    sage: T.group()
    Congruence Subgroup Gamma0(11)
    sage: T.dimension()
    2
    sage: T.cuspidal_subspace()
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
    Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: T.eisenstein_subspace()
    Eisenstein subspace of dimension 1 of Modular Forms space of dimension 2
    for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: M = ModularSymbols(11); M
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
    0 over Rational Field
    sage: M.weight()
    2
    sage: M.basis()
    ((1,0), (1,8), (1,9))
    sage: M.sign()
    0

设 `T_p` 表示通常的 Hecke 算子 (`p` 是质数)。
Hecke 算子 `T_2`, `T_3`, `T_5` 如何在模符号空间上作用？

.. link

::

    sage: M.T(2).matrix()
    [ 3  0 -1]
    [ 0 -2  0]
    [ 0  0 -2]
    sage: M.T(3).matrix()
    [ 4  0 -1]
    [ 0 -1  0]
    [ 0  0 -1]
    sage: M.T(5).matrix()
    [ 6  0 -1]
    [ 0  1  0]
    [ 0  0  1]

