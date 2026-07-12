********
接口问题
********

.. index::
   single: background, running Sage in

.. _section-background:

后台作业
========

没错，Sage 作业可以在 UNIX 系统上后台运行。通常做法是输入以下命令

.. code-block:: console

    $ nohup sage < command_file  > output_file &

使用 nohup 的优势在于即使你注销后，Sage 仍会继续运行。

目前在 (unix) ``top`` 命令的输出中，Sage 会显示为 "sage-ipython" 或 "python"，
但在 Sage 的未来版本中将显示为 ``sage``。

.. index::
   pair: referencing; Sage

引用 Sage
=========

请参阅 `citing Sage <https://doc.sagemath.org/html/en/faq/faq-general.html#i-want-to-cite-sage-in-a-publication-how-do-i-do-it>`_.

将 Sage 会话保存到日志
======================

没错，你可以将会话保存到日志。

(a) 你可以通过在后台运行 Sage (:ref:`section-background`)
将输出写入文件。

(b) 在 KDE konsole 中启动（仅适用于 Linux）。
进入 ``Settings`` `\rightarrow` ``History ...`` 并选择“无限制”。
启动你的会话。准备就绪后，进入 ``edit`` `\rightarrow` ``save history as ...``。

某些接口（如 Singular 或 GAP 的接口）允许你创建日志文件。
对于 Singular，有一个日志文件选项（在 ``singular.py`` 中）。
在 GAP 中，使用命令 ``LogTo``。

.. index:: LaTeX output

LaTeX 转换
==========

没错，你可以将某些结果输出为 LaTeX。

::

    sage: M = MatrixSpace(RealField(),3,3)
    sage: A = M([1,2,3, 4,5,6, 7,8,9])
    sage: print(latex(A))
    \left(\begin{array}{rrr}
        1.00000000000000 & 2.00000000000000 & 3.00000000000000 \\
        4.00000000000000 & 5.00000000000000 & 6.00000000000000 \\
        7.00000000000000 & 8.00000000000000 & 9.00000000000000
        \end{array}\right)

.. skip

::

    sage: view(A)

此时会自动调用 dvi 预览，在单独的窗口中显示生成的 LaTeX 输出。

还可以预览多元多项式和有理函数的 LaTeX：

::

    sage: x = PolynomialRing(QQ,3, 'x').gens()
    sage: f = x[0] + x[1] - 2*x[1]*x[2]
    sage: h = f /(x[1] + x[2])
    sage: print(latex(h))
    \frac{-2 x_{1} x_{2} + x_{0} + x_{1}}{x_{1} + x_{2}}

Sage 与其他计算机代数系统
=========================

如果 ``foo`` 是 Pari、GAP（不以分号结尾）、Singular、Maxima 命令，
对于 Pari，请输入 ``gp("foo")``，
对于 GAP，请输入 ``libgap.eval("foo")``，
对于 Singular，请输入 ``singular.eval("foo")``，
对于 Maxima，请输入 ``maxima("foo")``。
这些程序仅仅是将命令字符串发送到外部程序，并执行它，然后将结果读取回 Sage。
因此，如果外部程序没有安装并包含在 PATH 中，这些程序将无法工作。

.. index:: help in Sage

命令行上的 Sage 帮助
====================

如果你只知道 Sage 命令的部分名称，但想知道它在 Sage 中的位置，只需输入
``sage -grep <string>`` 即可在 Sage 源代码中查找所有出现 ``<string>`` 的地方。例如，

.. code-block:: console

    $ sage -grep berlekamp_massey
    matrix/all.py:from berlekamp_massey import berlekamp_massey
    matrix/berlekamp_massey.py:def berlekamp_massey(a):
    matrix/matrix.py:import berlekamp_massey
    matrix/matrix.py:            g =
    berlekamp_massey.berlekamp_massey(cols[i].list())

输入 ``help(foo)`` 或 ``foo??`` 获取帮助，
输入 ``foo.[tab]`` 来搜索 Sage 命令。输入 ``help()`` 获取 Python 命令的帮助。

例如

.. CODE-BLOCK:: python

    help(Matrix)

会在新屏幕中返回

.. skip

.. CODE-BLOCK:: text

    Help on cython_function_or_method in module sage.matrix.constructor:

    matrix(*args, **kwds)
        matrix(*args, **kwds)
        File: sage/matrix/constructor.pyx (starting at line 21)

            Create a matrix.

            This implements the ``matrix`` constructor::

                sage: matrix([[1,2],[3,4]])
                [1 2]
                [3 4]

            It also contains methods to create special types of matrices, see
            ``matrix.[tab]`` for more options. For example::
    --More--

输入 q 返回 Sage 屏幕。

.. index:: importing into Sage

读取和导入文件到 Sage
=====================

导入到 Sage 的文件必须以 `.py`` 结尾，例如 ``foo.py``，并且包含合法的 Python 语法。
前文 :ref:`section-permutation` 中的魔方群是一个简单示例。

另一种读取文件的方法是使用 ``load`` 或 ``attach`` 命令。
创建一个名为 ``example.sage`` 的文件（位于 Sage 的主目录中），内容如下：

.. skip

.. CODE-BLOCK:: python

    print("Hello World")
    print(2^3)

.. index:: load into Sage

使用 ``load`` 命令读取并执行 ``example.sage`` 文件：

.. skip

::

    sage: load("example.sage")
    Hello World
    8

.. index:: attach into Sage

你也可以将 Sage 文件 ``attach`` 到正在运行的会话中：

.. skip

::

    sage: attach("example.sage")
    Hello World
    8

现在，如果你更改 ``example.sage`` 并在 Sage 中输入空行，
那么 ``example.sage`` 的内容将自动重新加载到 Sage 中：

.. skip

::

    sage: !emacs example.sage&     #change 2^3 to 2^4
    sage:                          #hit return
    ***************************************************
                    Reloading 'example.sage'
    ***************************************************
    Hello World
    16

.. index:: Python and Sage

Sage 命令的 Python 语言程序代码
==============================================

假设你想知道 Sage 命令中用于计算置换群中心的 Python 程序是什么。
可以使用 Sage 的帮助界面查找文件名：

.. skip

::

    sage: PermutationGroup.center?
    Type:           instancemethod
    Base Class:     <class 'instancemethod'>
    String Form:    <unbound method PermutationGroup.center>
    Namespace:      Interactive
    File:           /home/wdj/sage/local/lib/python2.4/site-packages/sage/groups/permgroup.py
    Definition:     PermutationGroup.center(self)

现在你知道该命令位于 ``permgroup.py`` 文件中，并且知道该 Python 模块的目录。你可以使用编辑器来阅读源代码。

.. index:: special functions in Sage

Sage 中的“特殊函数”
=====================

Sage 有许多特殊函数（请参见参考手册 http://doc.sagemath.org/html/en/reference/functions/ ），
并且大多数可以进行符号操作。如果尚未实现，则其他符号包可能具有此功能。

通过 Maxima，可以进行一些符号操作：

::

    sage: maxima.eval("f:bessel_y (v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'
    sage: maxima.eval("diff (jacobi_sn (u, m), u)")
    'jacobi_cn(u,m)*jacobi_dn(u,m)'
    sage: jsn = lambda x: jacobi("sn",x,1)
    sage: P = plot(jsn,0,1, plot_points=20); Q = plot(lambda x:bessel_Y( 1, x), 1/2,1)
    sage: show(P)
    sage: show(Q)

除了 ``maxima`` 外，``pari`` 和 ``octave`` 也有特殊函数（实际上，Sage 封装了一些 ``pari`` 的特殊函数）。

下面是使用 Sage 接口（位于 sage/interfaces/octave.py）和 ``octave``
(https://www.gnu.org/software/octave/doc/latest) 的示例。

::

    sage: octave("atanh(1.1)")   ## optional - octave
    (1.52226,1.5708)

下面是使用 Sage 接口调用 ``pari`` 特殊函数的示例。

::

    sage: pari('2+I').besselk(3)
    0.0455907718407551 + 0.0289192946582081*I
    sage: pari('2').besselk(3)
    0.0615104584717420


Sage 是什么？
=============

Sage 是一个用来进行数论、代数和几何计算的框架，最初设计用于椭圆曲线和模形式的计算。
其长期目标是使其更广泛地应用于代数、几何和数论。它是开源的，并根据 GPL 条款免费提供。
可以从参考手册中的章节标题大致了解 Sage 涵盖的主题。

.. index::
   pair: Sage; history

Sage 的历史
-----------

Sage 由 William Stein 于 2004 年秋在哈佛大学创立，
0.1 版于 2005 年 1 月发布。该版本包括 Pari，但不包括 GAP 或 Singular。
0.2 版于 3 月发布，
0.3 版于 4 月发布，
0.4 版于 7 月发布。在此期间，Sage 增加了对 Cremona 数据库、多元多项式和大型有限域的支持。此外，还编写了更多的文档。
0.5 beta 版于 8 月发布，
0.6 beta 版于 9 月发布，
0.7 版于同月晚些时候发布。在此期间，增加了对向量空间、环、模符号和 Windows 用户的支持。
自 2005 年 10 月发布 0.8 版本以来，Sage 包含了 GAP 和 Singular 的完整发行版，尽管某些 GAP 数据库需要单独添加。
添加 Singular 并不容易，因为从源代码编译 Singular 非常困难。
0.9 版于 11 月发布。该版本经历了 34 次发布！
自 0.9.34 版（肯定在 0.10.0 版之前）以来，Maxima 和 clisp 被包含在 Sage 中。
0.10.0 版于 2006 年 1 月 12 日发布。
Sage 1.0 版于 2006 年 2 月初发布。
截至 2008 年 2 月，最新版本为 2.10.2。

许多人贡献了重要的代码和其他专业技术，例如协助在各种操作系统上进行编译。
一般来说，代码作者会在其 Python 文档的 AUTHOR 部分以及 Sage 网站的鸣谢部分予以体现。
