.. _chapter-plot:

****
绘图
****

Sage 可以使用 matplotlib, openmath, gnuplot 或 surf 进行绘图，
但在标准发行版中只有 matplotlib 和 openmath 包含在内。
有关 surf 的示例，请参见 :ref:`section-surface`。

在 Sage 中有许多不同的绘图方式。你可以使用 gnuplot 绘制函数（二维或三维）或者点集（仅支持二维），
你可以使用 Maxima 绘制微分方程的解（调用 gnuplot 或 openmath），
或者使用 Singular 的接口与绘图包 surf 进行绘图（不包含在 Sage 中）。
gnuplot 没有隐式绘图命令，所以如果你想使用隐式绘图绘制曲线或曲面，
如在章节‘ch:AG，代数几何’中所述，最好使用 Singular 的 surf 接口。


.. _section-piecewise:

二维函数绘图
============

默认的绘图方法使用优秀的 ``matplotlib`` 包。

要查看绘图，请输入 ``P.save("<path>/myplot.png")``，然后在诸如 gimp 之类的图形查看器中打开它。

如下例所示，可以绘制分段函数：

::

    sage: f1 = 1
    sage: f2 = 1-x
    sage: f3 = exp(x)
    sage: f4 = sin(2*x)
    sage: f = piecewise([((0,1),f1), ((1,2),f2), ((2,3),f3), ((3,10),f4)])
    sage: f.plot(x,0,10)
    Graphics object consisting of 1 graphics primitive

还可以生成其他函数图像：

雅可比椭圆函数
`\text{sn}(x,2)`, `-3<x<3` 的红色图像（请勿输入 ``....:``）:

::

    sage: L = [(i/100.0, maxima.eval('jacobi_sn (%s/100.0,2.0)'%i))
    ....:     for i in range(-300,300)]
    sage: show(line(L, rgbcolor=(3/4,1/4,1/8)))

`J`-贝塞尔函数 `J_2(x)`, `0<x<10` 的红色图像：

::

    sage: L = [(i/10.0, maxima.eval('bessel_j (2,%s/10.0)'%i)) for i in range(100)]
    sage: show(line(L, rgbcolor=(3/4,1/4,5/8)))

黎曼 `\zeta` 函数 `\zeta(1/2 + it)`, `0<t<30` 的紫色图像：

::

    sage: I = CDF.0
    sage: show(line([zeta(1/2 + k*I/6) for k in range(180)], rgbcolor=(3/4,1/2,5/8)))

.. _section-curve:

绘制曲线
========

在 Sage 中绘制曲线，可以使用 Singular 和 surf
(http://surf.sourceforge.net/ 也可以作为实验包使用)
或使用 matplotlib (Sage 自带)。

matplotlib
----------

这里有几个例子。要查看绘图，请输入
``p.save("<path>/my_plot.png")`` (其中 ``<path>`` 是你拥有写权限且希望保存绘图的目录路径)
并在查看器中查看 (例如 GIMP)。

蓝色 Nicomedes 共轭线：

::

    sage: L = [[1+5*cos(pi/2+pi*i/100), tan(pi/2+pi*i/100)*
    ....:     (1+5*cos(pi/2+pi*i/100))] for i in range(1,100)]
    sage: line(L, rgbcolor=(1/4,1/8,3/4))
    Graphics object consisting of 1 graphics primitive

蓝色（三叶）内旋轮线 (hypotrochoid)：

::

    sage: n = 4; h = 3; b = 2
    sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),
    ....:     n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
    sage: line(L, rgbcolor=(1/4,1/4,3/4))
    Graphics object consisting of 1 graphics primitive

蓝色（四叶）内旋轮线 (hypotrochoid)：

::

    sage: n = 6; h = 5; b = 2
    sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),
    ....:     n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
    sage: line(L, rgbcolor=(1/4,1/4,3/4))
    Graphics object consisting of 1 graphics primitive

红色帕斯卡蜗线 (limaçon of Pascal)：

::

    sage: L = [[sin(pi*i/100)+sin(pi*i/50),-(1+cos(pi*i/100)+cos(pi*i/50))]
    ....:     for i in range(-100,101)]
    sage: line(L, rgbcolor=(1,1/4,1/2))
    Graphics object consisting of 1 graphics primitive

浅绿色麦克劳林三等分角线 (trisectrix of Maclaurin)：

::

    sage: L = [[2*(1-4*cos(-pi/2+pi*i/100)^2),10*tan(-pi/2+pi*i/100)*
    ....:     (1-4*cos(-pi/2+pi*i/100)^2)] for i in range(1,100)]
    sage: line(L, rgbcolor=(1/4,1,1/8))
    Graphics object consisting of 1 graphics primitive


浅绿色伯努利双纽线 (lemniscate of Bernoulli，我们忽略 i=100，因为这样会导致除 0 错误)：

::

    sage: v = [(1/cos(-pi/2+pi*i/100), tan(-pi/2+pi*i/100)) for i in range(1,200) if i!=100 ]
    sage: L = [(a/(a^2+b^2), b/(a^2+b^2)) for a,b in v]
    sage: line(L, rgbcolor=(1/4,3/4,1/8))
    Graphics object consisting of 1 graphics primitive


.. index:: plot;curve using surf

surf
----

由于 ``surf`` 仅在类 UNIX 操作系统上可用（并且不包含在 Sage 中），
因此在 Sage 中使用以下命令绘图仅在类 UNIX 操作系统上可用。
顺带一提，surf 已包含在几个流行的 Linux 发行版中。


.. skip

::

    sage: s = singular.eval
    sage: s('LIB "surf.lib";')
    ...
    sage: s("ring rr0 = 0,(x1,x2),dp;")
    ''
    sage: s("ideal I = x1^3 - x2^2;")
    ''
    sage: s("plot(I);")
    ...

在 surf 窗口处于活动状态时按下 ``q`` 退出 surf 并返回 Sage。

你可以将此绘图保存为 surf 脚本。在弹出的 surf 窗口中，只需选择 ``file``, ``save as`` 等等。
（输入 ``q`` 或选择 ``file``, ``quit`` 关闭窗口。）

生成的图省略，但鼓励读者亲自尝试。

::

    s = singular
    s('LIB "surf.lib";')
    s("ring rr0 = 0,(x1,x2),dp;")
    s("ideal I = x13 - x22;")
    s("plot(I);")
    s('ring rr1 = 0,(x,y,z),dp;')
    s('ideal I(1) = 2x2-1/2x3 +1-y+1;')
    s('plot(I(1));')
    s('poly logo = ((x+3)3 + 2\*(x+3)2 - y2)\*(x3 -y2)\*((x-3)3-2\*(x-3)2-y2);')
    s('plot(logo);') Steiner surface
    s('ideal J(2) = x2\*y2+x2\*z2+y2\*z2-17\*x\*y\*z;')
    s('plot(J(2));')

openmath
========

Openmath 是由 W. Schelter 编写的 TCL/Tk GUI 绘图程序。

以下命令绘制函数
`\cos(2x)+2e^{-x}`

::

    sage: maxima.plot2d('cos(2*x) + 2*exp(-x)','[x,0,1]',  # not tested (pops up a window)
    ....:     '[plot_format,openmath]')

(Mac OS X 用户请注意：这些 ``openmath`` 命令是在一个 xterm shell 会话中运行的，而不是使用标准的 Mac 终端应用程序。)

::

    sage: maxima.eval('load("plotdf");')
    '".../share/maxima.../share/dynamics/plotdf.lisp"'
    sage: maxima.eval('plotdf(x+y,[trajectory_at,2,-0.1]); ')  # not tested

这里绘制了一个方向场（plotdf Maxima 包也由 W. Schelter 编写。）

多个函数的二维绘图：

::

    sage: maxima.plot2d('[x,x^2,x^3]','[x,-1,1]','[plot_format,openmath]')  # not tested

Openmath 还可以绘制 `z=f(x,y)` 形式的三维曲面，其中 `x` 和 `y` 在矩形范围内。
例如，这里有一个可以用鼠标移动的“实时”三维图：

::

    sage: maxima.plot3d ("sin(x^2 + y^2)", "[x, -3, 3]", "[y, -3, 3]",  # not tested
    ....:     '[plot_format, openmath]')

通过适当旋转，你可以查看轮廓图。

Tachyon 3D 绘图
===============

光线追踪包 Tachyon 随 Sage 一起分发。3D 图看起来非常漂亮，但通常需要更多的配置。
下面是一个参数空间曲线的示例：

::

    sage: f = lambda t: (t,t^2,t^3)
    sage: t = Tachyon(camera_center=(5,0,4))
    sage: t.texture('t')
    sage: t.light((-20,-20,40), 0.2, (1,1,1))
    sage: t.parametric_plot(f,-5,5,'t',min_depth=6)

输入 ``t.show()`` 查看该曲线。

参考手册中有其他示例。

gnuplot
=======

必须安装 ``gnuplot`` 才能运行这些命令。

.. index:: plot; a function

首先，这是绘制函数的方法：{plot!a function}

.. skip

::

    sage: maxima.plot2d('sin(x)','[x,-5,5]')
    sage: opts = '[gnuplot_term, ps], [gnuplot_out_file, "sin-plot.eps"]'
    sage: maxima.plot2d('sin(x)','[x,-5,5]',opts)
    sage: opts = '[gnuplot_term, ps], [gnuplot_out_file, "/tmp/sin-plot.eps"]'
    sage: maxima.plot2d('sin(x)','[x,-5,5]',opts)

默认情况下 eps 文件保存在当前目录，但你可以指定保存路径。


.. index:: plot; a parametric curve

下面是在平面上绘制参数曲线的示例：

.. skip

::

    sage: maxima.plot2d_parametric(["sin(t)","cos(t)"], "t",[-3.1,3.1])
    sage: opts = '[gnuplot_preamble, "set nokey"], [gnuplot_term, ps],
    ....:     [gnuplot_out_file, "circle-plot.eps"]'
    sage: maxima.plot2d_parametric(["sin(t)","cos(t)"], "t", [-3.1,3.1], options=opts)

下面是在三维空间中绘制参数曲面的示例:
{plot!a parametric surface}

.. skip

::

    sage: maxima.plot3d_parametric(["v*sin(u)","v*cos(u)","v"], ["u","v"],
    ....:     [-3.2,3.2],[0,3])
    sage: opts = '[gnuplot_term, ps], [gnuplot_out_file, "sin-cos-plot.eps"]'
    sage: maxima.plot3d_parametric(["v*sin(u)","v*cos(u)","v"], ["u","v"],
    ....:     [-3.2,3.2],[0,3],opts)

为了演示如何传递 gnuplot 选项，下面是关于绘制黎曼 `\zeta` 函数 `\zeta(s)` 上一系列点的示例
（使用 Pari 计算，但使用 Maxima 和 Gnuplot 绘图）: {plot!points} {Riemann zeta function}

.. skip

::

    sage: zeta_ptsx = [ (pari(1/2 + i*I/10).zeta().real()).precision(1)
    ....:     for i in range(70,150)]
    sage: zeta_ptsy = [ (pari(1/2 + i*I/10).zeta().imag()).precision(1)
    ....:     for i in range(70,150)]
    sage: maxima.plot_list(zeta_ptsx, zeta_ptsy)  # optional -- pops up a window.
    sage: opts='[gnuplot_preamble, "set nokey"], [gnuplot_term, ps],
    ....:     [gnuplot_out_file, "zeta.eps"]'
    sage: maxima.plot_list(zeta_ptsx, zeta_ptsy, opts) # optional -- pops up a window.

.. _section-surface:

绘制曲面
========

绘制曲面与绘制曲线并无太大区别，尽管语法稍有不同。特别是，你需要加载 ``surf``。
{plot!surface using surf}

.. skip

::

    sage: singular.eval('ring rr1 = 0,(x,y,z),dp;')
    ''
    sage: singular.eval('ideal I(1) = 2x2-1/2x3 +1-y+1;')
    ''
    sage: singular.eval('plot(I(1));')
    ...
