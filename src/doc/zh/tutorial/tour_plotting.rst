.. sage-doctest: needs sage.plot sage.symbolic

.. _section-plot:

绘图
========

Sage 可以生成二维和三维图形。

二维图形
---------------------

在二维中，Sage 可以绘制圆、线和多边形；在直角坐标系中绘制函数图形；
还可以绘制极坐标图、轮廓图和矢量场图。本文档展示了若干这些图形的例子。
有关使用 Sage 绘图的更多例子，请参见 :ref:`section-systems` 和 :ref:`section-maxima`，
以及 `Sage Constructions <http://doc.sagemath.org/html/en/constructions/>`_ 文档。

该命令生成一个位于原点的半径为 1 的黄色圆：

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0))
    Graphics object consisting of 1 graphics primitive

你还可以生成一个填充的圆：

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0), fill=True)
    Graphics object consisting of 1 graphics primitive

你还可以通过将圆赋值给变量来创建圆；这样做不会将圆绘制出来：

::

    sage: c = circle((0,0), 1, rgbcolor=(1,1,0))

要想绘制它，可以使用 ``c.show()`` 或 ``show(c)``，如下所示：

.. link

::

    sage: c.show()

或者，使用 ``c.save('filename.png')`` 将绘图保存到给定文件。

现在，这些“圆”看起来更像椭圆，因为坐标轴的比例不同。你可以这样修复这个问题：

.. link

::

    sage: c.show(aspect_ratio=1)

命令 ``show(c, aspect_ratio=1)`` 可以完成同样的事情，
或者你可以使用 ``c.save('filename.png', aspect_ratio=1)`` 保存图片。

绘制基本函数很容易：

::

    sage: plot(cos, (-5,5))
    Graphics object consisting of 1 graphics primitive

一旦你指定了变量名称，你还可以创建参数化图形：

::

    sage: x = var('x')
    sage: parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    Graphics object consisting of 1 graphics primitive

一定要注意，只有当原点在图形的视图范围内时，图形的轴才会相交，并且对于非常大的数值可能会使用科学计数法：

::

    sage: plot(x^2,(x,300,500))
    Graphics object consisting of 1 graphics primitive

你可以通过将多个图形相加来将他们组合在一起：

::

    sage: x = var('x')
    sage: p1 = parametric_plot((cos(x),sin(x)),(x,0,2*pi),rgbcolor=hue(0.2))
    sage: p2 = parametric_plot((cos(x),sin(x)^2),(x,0,2*pi),rgbcolor=hue(0.4))
    sage: p3 = parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    sage: show(p1+p2+p3, axes=false)

生成填充形状的一个好方法是生成点列表（示例中的 ``L``），
然后使用 ``polygon`` 命令绘制由这些点构成边界的形状。
例如，下面是一个绿色的三角形：

::

    sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),
    ....:     2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,3/4,1/2))
    sage: p
    Graphics object consisting of 1 graphics primitive

输入 ``show(p, axes=false)`` 来查看没有任何坐标轴的图形。

你可以向图形中添加文本：

::

    sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),
    ....:     6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))
    sage: t = text("hypotrochoid", (5,4), rgbcolor=(1,0,0))
    sage: show(p+t)

微积分老师经常在黑板上绘制以下图形：
arcsin 的多个周期：
即 :math:`y=\sin(x)` 对于 :math:`x` 在 :math:`-2\pi` 和 :math:`2\pi` 区间的图像，
围绕 45 度线翻转。以下 Sage 命令构造此图形：

::

    sage: v = [(sin(x),x) for x in srange(-2*float(pi),2*float(pi),0.1)]
    sage: line(v)
    Graphics object consisting of 1 graphics primitive

由于正切函数的值域比正弦函数大得多，如果你使用相同技巧绘制反正切的图像，你应该更改 *x* 轴的最大和最小坐标：

::

    sage: v = [(tan(x),x) for x in srange(-2*float(pi),2*float(pi),0.01)]
    sage: show(line(v), xmin=-20, xmax=20)

Sage 还能计算极坐标图、轮廓图和矢量场图（针对特殊类型的函数）。这里是一个轮廓图的例子：

::

    sage: f = lambda x,y: cos(x*y)
    sage: contour_plot(f, (-4, 4), (-4, 4))
    Graphics object consisting of 1 graphics primitive

三维图形
-----------------------

Sage 还可以用于创建三维图形。
在 notebook 和 REPL 中，这些图形将默认使用开源软件包 [ThreeJS]_ 显示，
该软件包支持使用鼠标交互式旋转和缩放图形。

使用 ``plot3d`` 绘制形如 `f(x, y) = z` 的函数图像：

::

    sage: x, y = var('x,y')
    sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
    Graphics3d Object

或者，你可以使用 ``parametric_plot3d`` 绘制参数曲面，
其中每个 `x, y, z` 由一个或两个变量（通常是 `u` 和 `v`）的函数确定。
前面的图形可以参数化地表达如下：

::

    sage: u, v = var('u, v')
    sage: f_x(u, v) = u
    sage: f_y(u, v) = v
    sage: f_z(u, v) = u^2 + v^2
    sage: parametric_plot3d([f_x, f_y, f_z], (u, -2, 2), (v, -2, 2))
    Graphics3d Object

在 Sage 中绘制 3D 曲面的第三种方法是 `implicit_plot3d``，
它绘制形如 `f(x, y, z) = 0` 的函数的轮廓（这定义了一组点）。
我们使用经典公式绘制一个球体：

::

    sage: x, y, z = var('x, y, z')
    sage: implicit_plot3d(x^2 + y^2 + z^2 - 4, (x,-2, 2), (y,-2, 2), (z,-2, 2))
    Graphics3d Object

下面是更多的例子：

`Yellow Whitney's umbrella <http://en.wikipedia.org/wiki/Whitney_umbrella>`__:

::

    sage: u, v = var('u,v')
    sage: fx = u*v
    sage: fy = u
    sage: fz = v^2
    sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -1, 1),
    ....:   frame=False, color="yellow")
    Graphics3d Object

`Cross cap <http://en.wikipedia.org/wiki/Cross-cap>`__:

::

    sage: u, v = var('u,v')
    sage: fx = (1+cos(v))*cos(u)
    sage: fy = (1+cos(v))*sin(u)
    sage: fz = -tanh((2/3)*(u-pi))*sin(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....:   frame=False, color="red")
    Graphics3d Object

挠环面：

::

    sage: u, v = var('u,v')
    sage: fx = (3+sin(v)+cos(u))*cos(2*v)
    sage: fy = (3+sin(v)+cos(u))*sin(2*v)
    sage: fz = sin(u)+2*cos(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....:   frame=False, color="red")
    Graphics3d Object

双纽线：

::

    sage: x, y, z = var('x,y,z')
    sage: f(x, y, z) = 4*x^2 * (x^2 + y^2 + z^2 + z) + y^2 * (y^2 + z^2 - 1)
    sage: implicit_plot3d(f, (x, -0.5, 0.5), (y, -1, 1), (z, -1, 1))
    Graphics3d Object
