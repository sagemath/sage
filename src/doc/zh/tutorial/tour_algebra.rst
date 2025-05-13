基本代数和微积分
==========================

Sage 能够进行多种与基本代数和微积分相关的计算，例如求解方程、微分、积分和拉普拉斯变换。
更多示例，请参阅
`Sage Constructions <http://doc.sagemath.org/html/en/constructions/>`_
。

在所有这些示例中，函数中的变量都需要使用 ``var(...)`` 定义。例如：

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

如果遇到 :class:`NameError` 错误，请检查是否拼写错误，或者是否忘记使用 ``var(...)`` 定义变量。


求解方程
-----------------

精确求解方程
~~~~~~~~~~~~~~~~~~~~~~~~~

``solve`` 函数用于求解方程。使用时，首先定义变量；
然后将方程（或方程组）和需要求解的变量作为 ``solve`` 的参数：

::

    sage: x = var('x')
    sage: solve(x^2 + 3*x + 2, x)
    [x == -2, x == -1]

你可以求解一元方程，其他变量作为参数：

::

    sage: x, b, c = var('x b c')
    sage: solve([x^2 + b*x + c == 0],x)
    [x == -1/2*b - 1/2*sqrt(b^2 - 4*c), x == -1/2*b + 1/2*sqrt(b^2 - 4*c)]

你也可以求解多元方程：

::

    sage: x, y = var('x, y')
    sage: solve([x+y==6, x-y==4], x, y)
    [[x == 5, y == 1]]

以下是由 Jason Grout 提供的使用 Sage 求解非线性方程组的示例：
首先，我们符号化地求解该方程组：

::

    sage: var('x y p q')
    (x, y, p, q)
    sage: eq1 = p+q==9
    sage: eq2 = q*y+p*x==-6
    sage: eq3 = q*y^2+p*x^2==24
    sage: solve([eq1,eq2,eq3,p==1],p,q,x,y)
    [[p == 1, q == 8, x == -4/3*sqrt(10) - 2/3, y == 1/6*sqrt(10) - 2/3], [p == 1, q == 8, x == 4/3*sqrt(10) - 2/3, y == -1/6*sqrt(10) - 2/3]]

对于解的数值近似，可以使用：

.. link

::

    sage: solns = solve([eq1,eq2,eq3,p==1],p,q,x,y, solution_dict=True)
    sage: [[s[p].n(30), s[q].n(30), s[x].n(30), s[y].n(30)] for s in solns]
    [[1.0000000, 8.0000000, -4.8830369, -0.13962039],
     [1.0000000, 8.0000000, 3.5497035, -1.1937129]]

（函数 ``n`` 用于打印数值近似，参数是精度的位数。）

数值求解方程
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

很多时候，``solve`` 无法找到指定方程或方程组的精确解。
此时可以使用 ``find_root`` 找到数值解。
例如，solve 对以下方程没有返回任何有意义的结果::

    sage: theta = var('theta')
    sage: solve(cos(theta)==sin(theta), theta)
    [sin(theta) == cos(theta)]

另一方面，可以使用 ``find_root`` 在区间 :math:`0 < \phi < \pi/2` 内找到上述方程的解::

    sage: phi = var('phi')
    sage: find_root(cos(phi)==sin(phi),0,pi/2)
    0.785398163397448...

微分、积分及其他
----------------------------------

Sage 可以对许多函数进行微分和积分。
例如，对 :math:`\sin(u)` 相对于 :math:`u` 进行微分，可以这样做：

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

计算 :math:`\sin(x^2)` 的四阶导数：

::

    sage: diff(sin(x^2), x, 4)
    16*x^4*sin(x^2) - 48*x^2*cos(x^2) - 12*sin(x^2)

分别计算 :math:`x^2+17y^2` 相对于 `x` 和 `y` 的偏导数：

::

    sage: x, y = var('x,y')
    sage: f = x^2 + 17*y^2
    sage: f.diff(x)
    2*x
    sage: f.diff(y)
    34*y

接下来讨论积分，包括不定积分和定积分。计算
:math:`\int x\sin(x^2)\, dx` 和
:math:`\int_0^1 \frac{x}{x^2+1}\, dx`

::

    sage: integral(x*sin(x^2), x)
    -1/2*cos(x^2)
    sage: integral(x/(x^2+1), x, 0, 1)
    1/2*log(2)

计算 :math:`\frac{1}{x^2-1}` 的部分分式分解：

::

    sage: f = 1/((1+x)*(x-1))
    sage: f.partial_fraction(x)
    -1/2/(x + 1) + 1/2/(x - 1)

.. _section-systems:

求解微分方程
------------------------------

你可以用 Sage 来研究常微分方程。
求解方程 :math:`x'+x-1=0`：

::

    sage: t = var('t')    # define a variable t
    sage: x = function('x')(t)   # define x to be a function of that variable
    sage: DE = diff(x, t) + x - 1
    sage: desolve(DE, [x,t])
    (_C + e^t)*e^(-t)

这里使用 Sage 与 Maxima [Max]_ 的接口，因此其输出可能与其他 Sage 输出有所不同。
上面示例中，输出表示该微分方程的一般解是
:math:`x(t) = e^{-t}(e^{t}+c)`。

你还可以计算拉普拉斯变换；
计算 :math:`t^2e^t -\sin(t)` 的拉普拉斯变换如下：

::

    sage: s = var("s")
    sage: t = var("t")
    sage: f = t^2*exp(t) - sin(t)
    sage: f.laplace(t,s)
    -1/(s^2 + 1) + 2/(s - 1)^3

这里是一个更复杂的示例。左侧连接到墙上的耦合弹簧的平衡位移

.. CODE-BLOCK:: text

    |------\/\/\/\/\---|mass1|----\/\/\/\/\/----|mass2|
             spring1               spring2

由二阶微分方程组建模

.. math::

    m_1 x_1'' + (k_1+k_2) x_1 - k_2 x_2 = 0

    m_2 x_2''+ k_2 (x_2-x_1) = 0,


其中 :math:`m_{i}` 是物体 *i* 的质量，:math:`x_{i}` 是质量 *i* 的平衡位移，:math:`k_{i}` 是弹簧 *i* 的弹簧常数。


**示例：** 使用 Sage 求解上述问题，其中
:math:`m_{1}=2`, :math:`m_{2}=1`, :math:`k_{1}=4`,
:math:`k_{2}=2`, :math:`x_{1}(0)=3`, :math:`x_{1}'(0)=0`,
:math:`x_{2}(0)=3`, :math:`x_{2}'(0)=0`.

解：对第一个方程进行拉普拉斯变换（符号 :math:`x=x_{1}`, :math:`y=x_{2}`）：

::

    sage: t,s = SR.var('t,s')
    sage: x = function('x')
    sage: y = function('y')
    sage: f = 2*x(t).diff(t,2) + 6*x(t) - 2*y(t)
    sage: f.laplace(t,s)
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

输出虽然难以阅读，但其表示

.. math:: -2x'(0) + 2s^2 \cdot X(s) - 2sx(0) - 2Y(s) + 6X(s) = 0


（其中小写函数如 :math:`x(t)` 的拉普拉斯变换是大写函数 :math:`X(s)`）。
对第二个方程进行拉普拉斯变换：

::

    sage: de2 = maxima("diff(y(t),t, 2) + 2*y(t) - 2*x(t)")
    sage: lde2 = de2.laplace("t","s"); lde2.sage()
    s^2*laplace(y(t), t, s) - s*y(0) - 2*laplace(x(t), t, s) + 2*laplace(y(t), t, s) - D[0](y)(0)

这表示

.. math:: -Y'(0) + s^2Y(s) + 2Y(s) - 2X(s) - sy(0) = 0.

代入初始条件 :math:`x(0)`, :math:`x'(0)`, :math:`y(0)`, 和 :math:`y'(0)`，
并求解所得的两个方程：

::

    sage: var('s X Y')
    (s, X, Y)
    sage: eqns = [(2*s^2+6)*X-2*Y == 6*s, -2*X +(s^2+2)*Y == 3*s]
    sage: solve(eqns, X,Y)
    [[X == 3*(s^3 + 3*s)/(s^4 + 5*s^2 + 4),
      Y == 3*(s^3 + 5*s)/(s^4 + 5*s^2 + 4)]]

此时进行逆拉普拉斯变换即可得到答案：

::

    sage: var('s t')
    (s, t)
    sage: inverse_laplace((3*s^3 + 9*s)/(s^4 + 5*s^2 + 4),s,t)
    cos(2*t) + 2*cos(t)
    sage: inverse_laplace((3*s^3 + 15*s)/(s^4 + 5*s^2 + 4),s,t)
    -cos(2*t) + 4*cos(t)

因此，解为

.. math:: x_1(t) = \cos(2t) + 2\cos(t), \quad x_2(t) = 4\cos(t) - \cos(2t).

可以使用参数方式绘制函数图像

::

    sage: t = var('t')
    sage: P = parametric_plot((cos(2*t) + 2*cos(t), 4*cos(t) - cos(2*t) ),
    ....:     (t, 0, 2*pi), rgbcolor=hue(0.9))
    sage: show(P)

也可以分开绘制两个函数的图像

::

    sage: t = var('t')
    sage: p1 = plot(cos(2*t) + 2*cos(t), (t,0, 2*pi), rgbcolor=hue(0.3))
    sage: p2 = plot(4*cos(t) - cos(2*t), (t,0, 2*pi), rgbcolor=hue(0.6))
    sage: show(p1 + p2)

有关绘图的更多信息，请参见 :ref:`section-plot`。
有关微分方程的更多信息，请参见 [NagleEtAl2004]_ 的第 5.5 节。


欧拉法求解微分方程组
----------------------------------------------------

在下一个示例中，我们将演示欧拉法求解一阶和二阶常微分方程。
首先回顾一下一阶方程的基本思想。给定初值问题的形式为

.. math::

    y'=f(x,y), \quad y(a)=c,

我们要找到解在 :math:`x=b` 处的近似值，其中 :math:`b>a`。

回顾导数的定义

.. math::  y'(x) \approx \frac{y(x+h)-y(x)}{h},


其中 :math:`h>0` 是一个给定且极小的数。
结合微分方程可以得到 :math:`f(x,y(x))\approx \frac{y(x+h)-y(x)}{h}`。
现在求解 :math:`y(x+h)`:

.. math::   y(x+h) \approx y(x) + h\cdot f(x,y(x)).


如果我们把 :math:`h \cdot f(x,y(x))` 称为“校正项”（因为没有更好的名称）,
把 :math:`y(x)` 称为“`y` 的旧值”，
把 :math:`y(x+h)` 称为“`y` 的新值”，
那么这个近似可以重新表示为

.. math::   y_{new} \approx y_{old} + h\cdot f(x,y_{old}).


如果我们将从 `a` 到 `b` 的区间分成 `n` 步，
使得 :math:`h=\frac{b-a}{n}`，那么我们可以在表中记录此方法的信息。

============== =======================   =====================
:math:`x`      :math:`y`                 :math:`h\cdot f(x,y)`
============== =======================   =====================
:math:`a`      :math:`c`                 :math:`h\cdot f(a,c)`
:math:`a+h`    :math:`c+h\cdot f(a,c)`         ...
:math:`a+2h`   ...
...
:math:`b=a+nh` ???                             ...
============== =======================   =====================


我们的目标是逐行填满表中的所有空白，直到到达 ??? 条目，这就是欧拉法对 :math:`y(b)` 的近似值。

求解微分方程组的思想与之类似。

**示例：** 数值近似 :math:`z(t)` 在 :math:`t=1` 处的值，使用欧拉法的 4 个步骤，
其中 :math:`z''+tz'+z=0`, :math:`z(0)=1`, :math:`z'(0)=0`。

我们必须将二阶常微分方程简化为两个一阶常微分方程组（使用 :math:`x=z`, :math:`y=z'`）并应用欧拉法：

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

因此，:math:`z(1)\approx 0.65`.

我们还可以绘制点 :math:`(x,y)` 以获得曲线的近似图。
函数 ``eulers_method_2x2_plot`` 将执行此操作；
为了使用它，我们需要定义函数 `f` 和 `g`，
它们接受一个带有三个坐标的参数：(`t`, `x`,`y`)。

::

    sage: f = lambda z: z[2]        # f(t,x,y) = y
    sage: g = lambda z: -sin(z[1])  # g(t,x,y) = -sin(x)
    sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)

此时，``P`` 存储了两个图： ``P[0]``, `x` 相对于 `t` 的图, 以及 ``P[1]``, `y` 相对于 `t` 的图。
我们可以通过如下代码绘制这两个图：

.. link

::

    sage: show(P[0] + P[1])

（有关绘图的更多信息，请参见 :ref:`section-plot`。）

特殊函数
-----------------

Sage 利用 PARI [GAP]_ 和 Maxima [Max]_ ,实现了多种正交多项式和特殊函数。
这些函数在 Sage 参考手册的相应部分（“正交多项式”和“特殊函数”）中有详细文档。

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

此时，Sage 仅将这些函数包装用于数值使用。
对于符号使用，请直接使用 Maxima 接口，如以下示例：

::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'


向量微积分
---------------

参见
`Vector Calculus Tutorial <http://doc.sagemath.org/html/en/thematic_tutorials/vector_calculus.html>`__.
