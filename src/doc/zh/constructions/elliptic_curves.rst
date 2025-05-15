.. index:: elliptic curves

***************
椭圆曲线
***************

导子
=========

如何在 Sage 中计算椭圆曲线（在 :math:`\QQ` 上）的导子？

在 Sage 中使用 ``EllipticCurve`` 命令定义椭圆曲线 :math:`E` 后，
导子就是与 :math:`E` 关联的若干“方法”之一。以下是语法示例（来自教程第 2.4 节“模形式”）：


::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E.conductor()
    10351

:math:`j`-不变量
=====================

如何在 Sage 中计算椭圆曲线的 :math:`j`-不变量？

与 ``EllipticCurve`` 类相关的其他方法包括 ``j_invariant``,
``discriminant`` 和 ``weierstrass_model``。以下是语法示例：

::

    sage: E = EllipticCurve([0, -1, 1, -10, -20])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
    sage: E.j_invariant()
    -122023936/161051
    sage: E.short_weierstrass_model()
    Elliptic Curve defined by y^2  = x^3 - 13392*x - 1080432 over Rational Field
    sage: E.discriminant()
    -161051
    sage: E = EllipticCurve(GF(5),[0, -1, 1, -10, -20])
    sage: E.short_weierstrass_model()
    Elliptic Curve defined by y^2  = x^3 + 3*x + 3 over Finite Field of size 5
    sage: E.j_invariant()
    4

.. index:: elliptic curves

E 上的 :math:`GF(q)`-有理点
========================================

如何计算有限域上椭圆曲线的点数？

给定一个定义在 :math:`\mathbb{F} = GF(q)` 上的椭圆曲线，
Sage 可以计算其 :math:`\mathbb{F}`-有理点集。

::

    sage: E = EllipticCurve(GF(5),[0, -1, 1, -10, -20])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + 4*x^2 over Finite Field of size 5
    sage: E.points()
    [(0 : 1 : 0), (0 : 0 : 1), (0 : 4 : 1), (1 : 0 : 1), (1 : 4 : 1)]
    sage: E.cardinality()
    5
    sage: G = E.abelian_group()
    sage: G
    Additive abelian group isomorphic to Z/5 embedded in Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 + 4*x^2 over Finite Field of size 5
    sage: G.permutation_group()
    Permutation Group with generators [(1,2,3,4,5)]

.. index::
   pair: modular form; elliptic curve

与 :math:`\QQ` 上椭圆曲线相关的模形式
========================================================================

设 :math:`E` 是一个“良好”的椭圆曲线，其方程具有整数系数。
设 :math:`N` 为 :math:`E` 的导子，并且对于每个 :math:`n`，
设 :math:`a_n` 是出现在 :math:`E` 的 Hasse-Weil :math:`L`-函数中的数字。
Taniyama-Shimura 猜想（已被 Wiles 证明）表明存在一个权重为 2、级别为 :math:`N` 的模形式，
它是 Hecke 算子下的特征形式，并具有傅里叶级数 :math:`\sum_{n = 0}^\infty a_n q^n`。
Sage 可以计算与 :math:`E` 相关的序列 :math:`a_n`。以下是一个示例。

::

    sage: E = EllipticCurve([0, -1, 1, -10, -20])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
    sage: E.conductor()
    11
    sage: E.anlist(20)
    [0, 1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0, 2]
    sage: E.analytic_rank()
    0
