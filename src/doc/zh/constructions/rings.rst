**
环
**

.. index::
   pair: matrix; ring

.. _section_matrix-ring:

矩阵环
======

如何在 Sage 中构建有限环上的矩阵环？``MatrixSpace`` 构造函数接受任意环作为基环。下面是语法示例：

::

    sage: R = IntegerModRing(51)
    sage: M = MatrixSpace(R,3,3)
    sage: M(0)
    [0 0 0]
    [0 0 0]
    [0 0 0]
    sage: M(1)
    [1 0 0]
    [0 1 0]
    [0 0 1]
    sage: 5*M(1)
    [5 0 0]
    [0 5 0]
    [0 0 5]

.. index::
   pair: polynomial; ring

.. _section-polynomial-ring:

多项式环
========

如何在 Sage 中构建有限域上的多项式环？
下面是一个示例：

::

    sage: R = PolynomialRing(GF(97),'x')
    sage: x = R.gen()
    sage: f = x^2+7
    sage: f in R
    True

下面是使用 Singular 接口的示例：

::

    sage: R = singular.ring(97, '(a,b,c,d)', 'lp')
    sage: I = singular.ideal(['a+b+c+d', 'ab+ad+bc+cd', 'abc+abd+acd+bcd', 'abcd-1'])
    sage: R
    polynomial ring, over a field, global ordering
    // coefficients: ZZ/97...
    // number of vars : 4
    //        block   1 : ordering lp
    //                  : names    a b c d
    //        block   2 : ordering C
    sage: I
    a+b+c+d,
    a*b+a*d+b*c+c*d,
    a*b*c+a*b*d+a*c*d+b*c*d,
    a*b*c*d-1

下面是另一个使用 GAP 的方法：

::

    sage: R = libgap.PolynomialRing(GF(97), 4); R
    GF(97)[x_1,x_2,x_3,x_4]
    sage: I = R.IndeterminatesOfPolynomialRing(); I
    [ x_1, x_2, x_3, x_4 ]
    sage: x1, x2, x3, x4 = I
    sage: f = x1*x2 + x3; f
    x_1*x_2+x_3
    sage: f.Value(I,[1,1,1,1])
    Z(97)^34

.. index:: p-adics

.. _section-padics:

`p`-进数
========

如何在 Sage 中构建 `p`-进数?
Sage 在这方面取得了很大的进展（参见 David Harvey 和 David Roe 在 SageDays 上的演讲）。
这里只给出一些简单的示例。

要计算 ``Qp`` 中整数的 ``Zp`` 环的特征和剩余类域，请使用以下示例所示语法。

::

    sage: K = Qp(3)
    sage: K.residue_class_field()
    Finite Field of size 3
    sage: K.residue_characteristic()
    3
    sage: a = K(1); a
    1 + O(3^20)
    sage: 82*a
    1 + 3^4 + O(3^20)
    sage: 12*a
    3 + 3^2 + O(3^21)
    sage: a in K
    True
    sage: b = 82*a
    sage: b^4
    1 + 3^4 + 3^5 + 2*3^9 + 3^12 + 3^13 + 3^16 + O(3^20)

.. index::
   pair: polynomial; quotient ring

多项式的商环
============

如何在 Sage 中构建商环？

我们创建商环 `GF(97)[x]/(x^3+7)`，并展示许多基本函数。

::

    sage: R = PolynomialRing(GF(97),'x')
    sage: x = R.gen()
    sage: S = R.quotient(x^3 + 7, 'a')
    sage: a = S.gen()
    sage: S
    Univariate Quotient Polynomial Ring in a over Finite Field of size 97 with
    modulus x^3 + 7
    sage: S.is_field()
    True
    sage: a in S
    True
    sage: x in S
    True
    sage: S.polynomial_ring()
    Univariate Polynomial Ring in x over Finite Field of size 97
    sage: S.modulus()
    x^3 + 7
    sage: S.degree()
    3

在 Sage 中，``in`` 表示存在对该环的“标准强制转换”。
因此整数 `x` 和 `a` 都在 `S` 中，
虽然 `x` 实际上需要被强制转换。

你还可以在商环中进行计算，而无需实际计算，然后使用 ``quo_rem`` 命令，如下所示：

::

    sage: R = PolynomialRing(GF(97),'x')
    sage: x = R.gen()
    sage: f = x^7+1
    sage: (f^3).quo_rem(x^7-1)
    (x^14 + 4*x^7 + 7, 8)
