******
微积分
******

这里有一些使用 Sage 进行微积分符号计算的例子。


.. index::
   pair: calculus; differentiation

微分
====

微分：

::

    sage: var('x k w')
    (x, k, w)
    sage: f = x^3 * e^(k*x) * sin(w*x); f
    x^3*e^(k*x)*sin(w*x)
    sage: f.diff(x)
    w*x^3*cos(w*x)*e^(k*x) + k*x^3*e^(k*x)*sin(w*x) + 3*x^2*e^(k*x)*sin(w*x)
    sage: latex(f.diff(x))
    w x^{3} \cos\left(w x\right) e^{\left(k x\right)} + k x^{3} e^{\left(k x\right)} \sin\left(w x\right) + 3 \, x^{2} e^{\left(k x\right)} \sin\left(w x\right)

如果你输入 ``view(f.diff(x))``，将会打开一个新的窗口显示编译后的输出。在 notebook 中，你可以在单元格中输入

.. CODE-BLOCK:: ipython

    var('x k w')
    f = x^3 * e^(k*x) * sin(w*x)
    show(f)
    show(f.diff(x))

然后按 ``shift-enter`` 来获得类似的结果。

你还可以在 notebook 单元格中，使用以下命令进行微分和积分

.. CODE-BLOCK:: ipython

    R = PolynomialRing(QQ,"x")
    x = R.gen()
    p = x^2 + 1
    show(p.derivative())
    show(p.integral())

或者在命令行中，使用以下命令进行微分和积分

::

    sage: R = PolynomialRing(QQ,"x")
    sage: x = R.gen()
    sage: p = x^2 + 1
    sage: p.derivative()
    2*x
    sage: p.integral()
    1/3*x^3 + x

此时你也可以输入 ``view(p.derivative())`` 或 ``view(p.integral())`` 以打开一个新窗口，
并由 LaTeX 对输出进行排版。

.. index::
   pair: calculus; critical points

临界点
------

你可以找到分段定义函数的临界点：

::

    sage: x = PolynomialRing(RationalField(), 'x').gen()
    sage: f1 = x^0
    sage: f2 = 1-x
    sage: f3 = 2*x
    sage: f4 = 10*x-x^2
    sage: f = piecewise([((0,1),f1), ((1,2),f2), ((2,3),f3), ((3,10),f4)])
    sage: f.critical_points()
    [5.0]

.. index:: Taylor series, power series

幂级数
------

Sage 提供了几种构建和处理幂级数的方法。

为了从函数表达式中获取泰勒级数，请在表达式上使用 ``.taylor()`` 方法::

    sage: var('f0 k x')
    (f0, k, x)
    sage: g = f0/sinh(k*x)^4
    sage: g.taylor(x, 0, 3)
    -62/945*f0*k^2*x^2 + 11/45*f0 - 2/3*f0/(k^2*x^2) + f0/(k^4*x^4)

可以通过 ``.series()`` 方法获得函数的形式幂级数展开::

    sage: (1/(2-cos(x))).series(x,7)
    1 + (-1/2)*x^2 + 7/24*x^4 + (-121/720)*x^6 + Order(x^7)

然而，目前对此类序列的某些操作很难执行。有两种替代方案：要么使用 Sage 的 Maxima 子系统以获取完整符号功能::

    sage: f = log(sin(x)/x)
    sage: f.taylor(x, 0, 10)
    -1/467775*x^10 - 1/37800*x^8 - 1/2835*x^6 - 1/180*x^4 - 1/6*x^2
    sage: maxima(f).powerseries(x,0)._sage_()
    sum(2^(2*i... - 1)*(-1)^i...*x^(2*i...)*bern(2*i...)/(i...*factorial(2*i...)), i..., 1, +Infinity)

要么你可以使用形式幂级数环进行快速计算。这些环缺乏符号函数::

    sage: R.<w> = QQ[[]]
    sage: ps = w + 17/2*w^2 + 15/4*w^4 + O(w^6); ps
    w + 17/2*w^2 + 15/4*w^4 + O(w^6)
    sage: ps.exp()
    1 + w + 9*w^2 + 26/3*w^3 + 265/6*w^4 + 413/10*w^5 + O(w^6)
    sage: (1+ps).log()
    w + 8*w^2 - 49/6*w^3 - 193/8*w^4 + 301/5*w^5 + O(w^6)
    sage: (ps^1000).coefficients()
    [1, 8500, 36088875, 102047312625, 1729600092867375/8]

.. index::
   pair: calculus; integration

积分
====

下面的 :ref:`section-riemannsums` 讨论了数值积分。

Sage 可以自行对一些简单函数进行积分：

::

    sage: f = x^3
    sage: f.integral(x)
    1/4*x^4
    sage: integral(x^3,x)
    1/4*x^4
    sage: f = x*sin(x^2)
    sage: integral(f,x)
    -1/2*cos(x^2)

Sage 还可以计算涉及极限的符号定积分。

::

    sage: var('x, k, w')
    (x, k, w)
    sage: f = x^3 * e^(k*x) * sin(w*x)
    sage: f.integrate(x)
    ((24*k^3*w - 24*k*w^3 - (k^6*w + 3*k^4*w^3 + 3*k^2*w^5 + w^7)*x^3 + 6*(k^5*w + 2*k^3*w^3 + k*w^5)*x^2 - 6*(3*k^4*w + 2*k^2*w^3 - w^5)*x)*cos(w*x)*e^(k*x) - (6*k^4 - 36*k^2*w^2 + 6*w^4 - (k^7 + 3*k^5*w^2 + 3*k^3*w^4 + k*w^6)*x^3 + 3*(k^6 + k^4*w^2 - k^2*w^4 - w^6)*x^2 - 6*(k^5 - 2*k^3*w^2 - 3*k*w^4)*x)*e^(k*x)*sin(w*x))/(k^8 + 4*k^6*w^2 + 6*k^4*w^4 + 4*k^2*w^6 + w^8)
    sage: integrate(1/x^2, x, 1, infinity)
    1


.. index: convolution

卷积
----

你可以计算任意分段函数与另一个函数的卷积（在定义域之外，它们被假定为零）。
以下是 `f`, `f*f` 和 `f*f*f` 的定义，
其中 `f(x)=1`, `0<x<1`:

::

    sage: x = PolynomialRing(QQ, 'x').gen()
    sage: f = piecewise([((0,1),1*x^0)])
    sage: g = f.convolution(f)
    sage: h = f.convolution(g)
    sage: set_verbose(-1)
    sage: P = f.plot(); Q = g.plot(rgbcolor=(1,1,0)); R = h.plot(rgbcolor=(0,1,1))

要查看此内容，请输入 ``show(P+Q+R)``。


.. _section-riemannsums:

黎曼和与梯形法积分
------------------

关于 `\int_a^bf(x)\, dx` 的数值近似，
其中 `f`  是分段函数，可以：


-  计算（用于绘图目的）根据梯形法则定义的分段线性函数，基于将其分割为 `N` 个子区间进行数值积分；

-  梯形法则给出的近似值；

-  计算（用于绘图目的）根据黎曼和（左端点、右端点或中点）定义的分段常数函数，
   基于将其分割为 `N` 个子区间进行数值积分；

-  黎曼和近似值给出的近似值。


::

    sage: f1(x) = x^2
    sage: f2(x) = 5-x^2
    sage: f = piecewise([[[0,1], f1], [RealSet.open_closed(1,2), f2]])
    sage: t = f.trapezoid(2); t
    piecewise(x|-->1/2*x on (0, 1/2), x|-->3/2*x - 1/2 on (1/2, 1), x|-->7/2*x - 5/2 on (1, 3/2), x|-->-7/2*x + 8 on (3/2, 2); x)
    sage: t.integral()
    piecewise(x|-->1/4*x^2 on (0, 1/2), x|-->3/4*x^2 - 1/2*x + 1/8 on (1/2, 1), x|-->7/4*x^2 - 5/2*x + 9/8 on (1, 3/2), x|-->-7/4*x^2 + 8*x - 27/4 on (3/2, 2); x)
    sage: t.integral(definite=True)
    9/4

.. index: Laplace transform

拉普拉斯变换
------------

如果你有一个分段定义的多项式函数，那么有一个“原生”命令用于计算拉普拉斯变换。
这将调用 Maxima，但值得注意的是，Maxima 无法（使用最后几个示例中的直接接口）处理这种类型的计算。

::

    sage: var('x s')
    (x, s)
    sage: f1(x) = 1
    sage: f2(x) = 1-x
    sage: f = piecewise([((0,1),f1), ((1,2),f2)])
    sage: f.laplace(x, s)
    -e^(-s)/s + (s + 1)*e^(-2*s)/s^2 + 1/s - e^(-s)/s^2

对于其他“合理”的函数，可以使用 Maxima 接口计算拉普拉斯变换：

::

    sage: var('k, s, t')
    (k, s, t)
    sage: f = 1/exp(k*t)
    sage: f.laplace(t,s)
    1/(k + s)

上面是计算拉普拉斯变换的一种方法

::

    sage: var('s, t')
    (s, t)
    sage: f = t^5*exp(t)*sin(t)
    sage: L = laplace(f, t, s); L
    3840*(s - 1)^5/(s^2 - 2*s + 2)^6 - 3840*(s - 1)^3/(s^2 - 2*s + 2)^5 +
    720*(s - 1)/(s^2 - 2*s + 2)^4

上面是另一种方法。

.. index:
   pair: differential equations; solve

常微分方程
==========

使用 Sage 接口与 Maxima 可以符号化地求解常微分方程。参见

.. skip

::

    sage: desolvers?

获取可用命令。
可以使用 Sage 接口与 Octave（一个实验性包）或 GSL（Gnu 科学库）中的例程来数值求解常微分方程。

例如，通过 Sage 的 Maxima 接口符号化地求解常微分方程（请勿输入 ``....:``）：

::

    sage: y=function('y')(x); desolve(diff(y,x,2) + 3*x == y, dvar = y, ics = [1,1,1])
    3*x - 2*e^(x - 1)
    sage: desolve(diff(y,x,2) + 3*x == y, dvar = y)
    _K2*e^(-x) + _K1*e^x + 3*x
    sage: desolve(diff(y,x) + 3*x == y, dvar = y)
    (3*(x + 1)*e^(-x) + _C)*e^x
    sage: desolve(diff(y,x) + 3*x == y, dvar = y, ics = [1,1]).expand()
    3*x - 5*e^(x - 1) + 3

    sage: f=function('f')(x); desolve_laplace(diff(f,x,2) == 2*diff(f,x)-f, dvar = f, ics = [0,1,2])
    x*e^x + e^x

    sage: desolve_laplace(diff(f,x,2) == 2*diff(f,x)-f, dvar = f)
    -x*e^x*f(0) + x*e^x*D[0](f)(0) + e^x*f(0)

.. index:
   pair: differential equations; plot

如果你已经安装了 ``Octave`` 和 ``gnuplot``，

::

    sage: octave.de_system_plot(['x+y','x-y'], [1,-1], [0,2]) # optional - octave

将在同一个图中绘制常微分方程组的两个图像 `(t,x(t)), (t,y(t))` （`t`-轴为横轴）

.. MATH::

    x' = x+y, x(0) = 1; y' = x-y, y(0) = -1,

对于 `0 \leq t \leq 2`。使用 ``desolve_system_rk4`` 也可以获得相同的结果::

    sage: x, y, t = var('x y t')
    sage: P=desolve_system_rk4([x+y, x-y], [x,y], ics=[0,1,-1], ivar=t, end_points=2)
    sage: p1 = list_plot([[i,j] for i,j,k in P], plotjoined=True)
    sage: p2 = list_plot([[i,k] for i,j,k in P], plotjoined=True, color='red')
    sage: p1+p2
    Graphics object consisting of 2 graphics primitives

该方程组也可以通过使用命令 ``desolve_system`` 来求解。

.. skip

::

    sage: t=var('t'); x=function('x',t); y=function('y',t)
    sage: des = [diff(x,t) == x+y, diff(y,t) == x-y]
    sage: desolve_system(des, [x,y], ics = [0, 1, -1])
    [x(t) == cosh(sqrt(2)*t), y(t) == sqrt(2)*sinh(sqrt(2)*t) - cosh(sqrt(2)*t)]

此命令的输出 *不* 是一对函数。

最后，可以使用幂级数求解线性微分方程：

::

    sage: R.<t> = PowerSeriesRing(QQ, default_prec=10)
    sage: a = 2 - 3*t + 4*t^2 + O(t^10)
    sage: b = 3 - 4*t^2 + O(t^7)
    sage: f = a.solve_linear_de(prec=5, b=b, f0=3/5)
    sage: f
    3/5 + 21/5*t + 33/10*t^2 - 38/15*t^3 + 11/24*t^4 + O(t^5)
    sage: f.derivative() - a*f - b
    O(t^4)

周期函数的傅里叶级数
====================

设 `f` 是一个周期为 `2L` 的实周期函数。
`f` 的傅里叶级数是

.. MATH::

   S(x) = \frac{a_0}{2} + \sum_{n=1}^\infty \left[a_n\cos\left(\frac{n\pi x}{L}\right) +
   b_n\sin\left(\frac{n\pi x}{L}\right)\right]

其中

.. MATH::

    a_n = \frac{1}{L}\int_{-L}^L
            f(x)\cos\left(\frac{n\pi x}{L}\right) dx,

并且

.. MATH::

    b_n = \frac{1}{L}\int_{-L}^L
            f(x)\sin\left(\frac{n\pi x}{L}\right) dx,

傅里叶系数 `a_n` 和 `b_n` 是通过声明 `f` 在一个周期内分段定义的函数并调用方法
``fourier_series_cosine_coefficient`` 和 ``fourier_series_sine_coefficient`` 来计算的，
而部分和是通过 ``fourier_series_partial_sum`` 获得的::

    sage: f = piecewise([((0,pi/2), -1), ((pi/2,pi), 2)])
    sage: f.fourier_series_cosine_coefficient(0)
    1
    sage: f.fourier_series_sine_coefficient(5)
    -6/5/pi
    sage: s5 = f.fourier_series_partial_sum(5); s5
    -6/5*sin(10*x)/pi - 2*sin(6*x)/pi - 6*sin(2*x)/pi + 1/2
    sage: plot(f, (0,pi)) + plot(s5, (x,0,pi), color='red')
    Graphics object consisting of 2 graphics primitives

.. PLOT::

    f = piecewise([((0,pi/2), -1), ((pi/2,pi), 2)])
    s5 = f.fourier_series_partial_sum(5)
    g = plot(f, (0,pi)) + plot(s5, (x,0,pi), color='red')
    sphinx_plot(g)
