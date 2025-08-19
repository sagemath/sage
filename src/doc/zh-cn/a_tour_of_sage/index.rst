.. _a-tour-of-sage:

===============
欢迎来到 Sage
===============

这是 Sage 作为计算器的简短介绍。

Sage 命令行的提示符为“``sage:``”。要尝试以下示例，只需输入提示符后的部分。

::

    sage: 3 + 5
    8

如果你在 Jupyter 笔记本上使用 Sage，那么同样地，将提示符后的所有内容放入输入单元格中，然后按 :kbd:`Shift-Enter` 键获得相应的输出。

插入符号 “^” 表示“乘方。

::

    sage: 57.1^100
    4.60904368661396e175

我们用 Sage 计算 :math:`2 \times 2` 矩阵的逆。

::

    sage: matrix([[1, 2], [3, 4]])^(-1)
    [  -2    1]
    [ 3/2 -1/2]

在这里，我们整合了一个简单的函数。

::

    sage: x = var('x')   # create a symbolic variable
    sage: integrate(sqrt(x) * sqrt(1 + x), x)
    1/4*((x + 1)^(3/2)/x^(3/2) + sqrt(x + 1)/sqrt(x))/((x + 1)^2/x^2 - 2*(x + 1)/x + 1)
    - 1/8*log(sqrt(x + 1)/sqrt(x) + 1) + 1/8*log(sqrt(x + 1)/sqrt(x) - 1)

这要求 Sage 解一元二次方程。符号 ``==`` 在 Sage 中代表相等。

::

    sage: a = var('a')
    sage: S = solve(x^2 + x == a, x); S
    [x == -1/2*sqrt(4*a + 1) - 1/2, x == 1/2*sqrt(4*a + 1) - 1/2]

结果是一份等式列表。

.. link

::

    sage: S[0].rhs()  # right hand side of the equation
    -1/2*sqrt(4*a + 1) - 1/2

当然，Sage 可以绘制各种有用的函数。

::

    sage: show(plot(sin(x) + sin(1.6*x), 0, 40))

.. image:: sin_plot.*


Sage 是一款非常强大的计算器。要体验它，首先我们创建一个 :math:`500 \times 500` 随机数矩阵。

::

    sage: m = random_matrix(RDF, 500)

用 Sage 计算矩阵的特征值并绘制它们只需要一秒钟的时间。

.. link

::

    sage: e = m.eigenvalues()  # about 1 second
    sage: w = [(i, abs(e[i])) for i in range(len(e))]
    sage: show(points(w))

.. image:: eigen_plot.*


Sage 可以处理非常大的数字，甚至是数百万或数十亿位的数字。

::

    sage: factorial(100)
    93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000

::

    sage: n = factorial(1000000)  # about 1 second
    sage: len(n.digits())
    5565709

这将计算 :math:`\pi` 的至少 100 位数字。

::

    sage: N(pi, digits=100)
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068

这要求 Sage 对两个变量的多项式进行因式分解。

::

    sage: R.<x,y> = QQ[]
    sage: F = factor(x^99 + y^99)
    sage: F
    (x + y) * (x^2 - x*y + y^2) * (x^6 - x^3*y^3 + y^6) *
    (x^10 - x^9*y + x^8*y^2 - x^7*y^3 + x^6*y^4 - x^5*y^5 +
     x^4*y^6 - x^3*y^7 + x^2*y^8 - x*y^9 + y^10) *
    (x^20 + x^19*y - x^17*y^3 - x^16*y^4 + x^14*y^6 + x^13*y^7 -
     x^11*y^9 - x^10*y^10 - x^9*y^11 + x^7*y^13 + x^6*y^14 -
     x^4*y^16 - x^3*y^17 + x*y^19 + y^20) * (x^60 + x^57*y^3 -
     x^51*y^9 - x^48*y^12 + x^42*y^18 + x^39*y^21 - x^33*y^27 -
     x^30*y^30 - x^27*y^33 + x^21*y^39 + x^18*y^42 - x^12*y^48 -
     x^9*y^51 + x^3*y^57 + y^60)
    sage: F.expand()
    x^99 + y^99

Sage 只需不到 1 秒钟的时间就能计算出将一亿分割为正整数之和的方法数。

::

    sage: z = Partitions(10^8).cardinality()  # about .1 second
    sage: z
    1760517045946249141360373894679135204009...

Sage 是世界上最先进的开源数学软件。
