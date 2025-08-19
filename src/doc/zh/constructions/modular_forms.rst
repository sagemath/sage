.. index:: modular forms

******
模形式
******

SageMath 的计算专长之一是模形式（非常技术性的领域），它能做的远远超出这篇非常简短的介绍所给出的内容。

尖点形式
========

如何使用 Sage 计算尖点形式空间的维数？

要计算 Gamma 尖点形式空间的维数，请使用命令 ``dimension_cusp_forms``。
以下是教程中“模形式”章节的一个示例：

::

    sage: from sage.modular.dims import dimension_cusp_forms
    sage: dimension_cusp_forms(Gamma0(11),2)
    1
    sage: dimension_cusp_forms(Gamma0(1),12)
    1
    sage: dimension_cusp_forms(Gamma1(389),2)
    6112

相关命令：``dimension_new__cusp_forms_gamma0`` （用于新形式的维数）、
``dimension_new__cusp_forms_gamma0`` （用于模形式）
以及 ``dimension_eis`` （用于艾森斯坦级数）。
这些命令的语法类似 - 可以参阅参考手册中的示例。

.. index:: cosets of Gamma_0

陪集表示
=====================

算术商 `H/\Gamma` 的基本定义域的显式表示可以通过 `\Gamma` 在 `SL_2(\ZZ)` 上的陪集来确定。
这些陪集是如何在 Sage 中计算的呢？

以下是计算 `SL_2(\ZZ)/\Gamma_0(11)` 的陪集表示的示例：

::

    sage: G = Gamma0(11); G
    Congruence Subgroup Gamma0(11)
    sage: list(G.coset_reps())
    [
    [1 0]  [ 0 -1]  [1 0]  [ 0 -1]  [ 0 -1]  [ 0 -1]  [ 0 -1]  [ 0 -1]
    [0 1], [ 1  0], [1 1], [ 1  2], [ 1  3], [ 1  4], [ 1  5], [ 1  6],
    <BLANKLINE>
    [ 0 -1]  [ 0 -1]  [ 0 -1]  [ 0 -1]
    [ 1  7], [ 1  8], [ 1  9], [ 1 10]
    ]


.. index:: modular symbols, Hecke operators

模符号和 Hecke 算子
===================

接下来我们展示如何在级别为 1、权重为 12 的模符号空间上计算 Hecke 算子。

::

    sage: M = ModularSymbols(1,12)
    sage: M.basis()
    ([X^8*Y^2,(0,0)], [X^9*Y,(0,0)], [X^10,(0,0)])
    sage: t2 = M.T(2)
    sage: f = t2.charpoly('x'); f
    x^3 - 2001*x^2 - 97776*x - 1180224
    sage: factor(f)
    (x - 2049) * (x + 24)^2
    sage: M.T(11).charpoly('x').factor()
    (x - 285311670612) * (x - 534612)^2

这里 ``t2`` 表示 `\QQ` 上对权重为 `12` 、符号为 `0`
且维度为 `3` 的 `\Gamma_0(1)` 的完整模符号空间的 Hecke 算子 `T_2`。

::

    sage: M = ModularSymbols(Gamma1(6),3,sign=0)
    sage: M
    Modular Symbols space of dimension 4 for Gamma_1(6) of weight 3 with sign 0
    over Rational Field
    sage: M.basis()
    ([X,(0,5)], [X,(3,5)], [X,(4,5)], [X,(5,5)])
    sage: M._compute_hecke_matrix_prime(2).charpoly()
    x^4 - 17*x^2 + 16
    sage: M.integral_structure()
    Free module of degree 4 and rank 4 over Integer Ring
    Echelon basis matrix:
    [1 0 0 0]
    [0 1 0 0]
    [0 0 1 0]
    [0 0 0 1]

更多示例请参阅教程或参考手册中的模形式章节。

亏格公式 (Genus formulas)
=========================

Sage 可以计算 `X_0(N)`, `X_1(N)` 及相关曲线的亏格 (genus)。
以下是一些语法示例：

::

    sage: from sage.modular.dims import dimension_cusp_forms
    sage: dimension_cusp_forms(Gamma0(22))
    2
    sage: dimension_cusp_forms(Gamma0(30))
    3
    sage: dimension_cusp_forms(Gamma1(30))
    9

请参阅计算模形式空间维数的代码（在 ``sage/modular/dims.py`` 中）
或 Oesterlé 和 Cohen 的论文 {CO}，获取详细信息。
