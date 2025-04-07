.. _chapter-interactive_shell:

*********************
交互式 Shell
*********************
在本教程的大部分内容中，我们假定你使用 ``sage`` 命令启动 Sage 解释器。
这将启动一个定制版的 IPython Shell，并导入许多函数和类，使它们可以直接从命令提示符使用。
可以通过编辑 ``$SAGE_ROOT/ipythonrc`` 文件进行进一步的自定义。
启动 Sage 后，会输出以下类似内容：

.. CODE-BLOCK:: text

    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘


    sage:

要退出 Sage 只需按 Ctrl-D 或输入 ``quit`` 或 ``exit``。

.. skip

::

    sage: quit
    Exiting Sage (CPU time 0m0.00s, Wall time 0m0.89s)

Wall time 指的是墙上的挂钟走过的时间。因为 CPU 时间不会跟踪子进程（如 GAP 或 Singular）消耗的时间。

（请避免在终端中使用 ``kill -9`` 杀死 Sage 进程，
因为 Sage 可能无法终止子进程，例如 Maple 进程，或清理 ``$HOME/.sage/tmp`` 中的临时文件。）

Sage 会话
=================

会话是从 Sage 启动到退出期间的输入输出序列。Sage 通过IPython 记录所有 Sage 输入。
实际上，如果你使用的是交互式 Shell（而不是 notebook 界面），
你可以随时输入 ``%history`` （或 ``%hist``）来列出迄今为止输入的所有命令行。
在 Sage 提示符下输入 ``?`` 可以了解有关 IPython 的更多信息，例如，
“IPython 提供带编号的提示符...并缓存输入和输出。所有输入都会保存，
并且可以作为变量检索（除了常用的箭头键召回外）。以下全局变量始终存在（所以不要覆盖它们！）”：

.. CODE-BLOCK:: text

      _:  上一次输入 (交互式 SHell 和 notebook 均适用)
      __: 上两次输入 (仅交互式 Shell 适用)
      _oh : 所有输入的列表 (仅交互式 Shell 适用)

例如：

.. skip

::

    sage: factor(100)
     _1 = 2^2 * 5^2
    sage: kronecker_symbol(3,5)
     _2 = -1
    sage: %hist   # This only works from the interactive shell, not the notebook.
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    sage: _oh
     _4 = {1: 2^2 * 5^2, 2: -1}
    sage: _i1
     _5 = 'factor(ZZ(100))\n'
    sage: eval(_i1)
     _6 = 2^2 * 5^2
    sage: %hist
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    4: _oh
    5: _i1
    6: eval(_i1)
    7: %hist

我们在本教程和其他 Sage 文档中均省略了输出编号。

你还可以在会话中将输入列表储存在宏中。

.. skip

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: M = ModularSymbols(37)
    sage: %hist
    1: E = EllipticCurve([1,2,3,4,5])
    2: M = ModularSymbols(37)
    3: %hist
    sage: %macro em 1-2
    Macro `em` created. To execute, type its name (without quotes).


.. skip

::

    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E = 5
    sage: M = None
    sage: em
    Executing Macro...
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field

在使用交互式 Shell 时，任何 UNIX Shell 命令都可以通过在 Sage 前面加上感叹号 ``!`` 来执行。例如：

.. skip

::

    sage: !ls
    auto  example.sage glossary.tex  t  tmp  tut.log  tut.tex

返回当前目录的列表。

``PATH`` 变量将 Sage 的 bin 目录放在最前端，
因此如果运行 ``gp``, ``gap``, ``singular``, ``maxima`` 等等，你会得到随 Sage 附带的版本。

.. skip

::

    sage: !gp
    Reading GPRC: /etc/gprc ...Done.

                               GP/PARI CALCULATOR Version 2.2.11 (alpha)
                      i686 running linux (ix86/GMP-4.1.4 kernel) 32-bit version
    ...
    sage: !singular
                         SINGULAR                             /  Development
     A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \

记录输入和输出
========================

记录 Sage 会话不同于保存会话（参见 :ref:`section-save`）。
要记录输入（和可选输出），请使用 ``logstart`` 命令。输入 ``logstart?`` 了解更多详情。
你可以使用这个命令记录你输入的所有内容、所有输出，甚至可以在未来的会话中重现输入（通过重新加载日志文件）。

.. skip

.. CODE-BLOCK:: shell-session

    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: logstart setup
    Activating auto-logging. Current session state plus future input saved.
    Filename       : setup
    Mode           : backup
    Output logging : False
    Timestamping   : False
    State          : active
    sage: E = EllipticCurve([1,2,3,4,5]).minimal_model()
    sage: F = QQ^3
    sage: x,y = QQ['x,y'].gens()
    sage: G = E.gens()
    sage:
    Exiting Sage (CPU time 0m0.61s, Wall time 0m50.39s).
    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: load("setup")
    Loading log file <setup> one line at a time...
    Finished replaying log file <setup>
    sage: E
    Elliptic Curve defined by y^2 + x*y  = x^3 - x^2 + 4*x + 3 over Rational
    Field
    sage: x*y
    x*y
    sage: G
    [(2 : 3 : 1)]

如果你在 Linux KDE 终端 ``konsole`` 中使用 Sage，那么可以按照以下步骤保存会话：
在 ``konsole`` 中启动 Sage 后，选择“设置”，然后“历史记录...”，然后“设置为无限制”。
当你准备保存会话时，选择“编辑”，然后“保存历史记录为...”，并输入一个名称将会话的文本保存到你的计算机。
保存这个文件后，你可以将其加载到编辑器（例如 xemacs）并打印出来。

粘贴忽略提示符
=====================

假设你正在阅读 Sage 或 Python 计算的会话，并希望将它们复制到 Sage 中。
但是有 ``>>>`` 或 ``sage:`` 提示符很烦人。实际上，你可以将包含提示符的示例复制并粘贴到 Sage 中。
换句话说，默认情况下，Sage 解析器在传递给 Python 之前会删除任何前导 ``>>>`` 或 ``sage:`` 提示符。例如：

.. skip

::

    sage: 2^10
    1024
    sage: sage: sage: 2^10
    1024
    sage: >>> 2^10
    1024

命令计时
===============

如果你在输入的开头放置 ``%time`` 命令，那么命令执行的时间将显示在输出后。
例如，我们可以比较几种幂运算的运行时间。
这些计时在你电脑上可能会有很大不同，甚至在不同版本的 Sage 之间也会有所不同。
首先是原生 Python：

.. skip

::

    sage: %time a = int(1938)^int(99484)
    CPU times: user 0.66 s, sys: 0.00 s, total: 0.66 s
    Wall time: 0.66

这意味着总共耗时 0.66 秒，"Wall time" 即墙上挂钟的时间为 0.66 秒。
如果你的计算机负载较重，wall time 可能比 CPU 时间长很多。

还可以使用 ``timeit`` 函数来尝试在大量迭代命令下获取时间。
这提供了稍微不同的信息，并且需要输入命令字符串来计时。

.. skip

::

    sage: timeit("int(1938)^int(99484)")
    5 loops, best of 3: 44.8 ms per loop

接下来我们使用原生 Sage Integer 类型，它是用 Cython 调用 GMP 库实现的：

.. skip

::

    sage: %time a = 1938^99484
    CPU times: user 0.04 s, sys: 0.00 s, total: 0.04 s
    Wall time: 0.04

使用 PARI 的 C 语言接口：

.. skip

::

    sage: %time a = pari(1938)^pari(99484)
    CPU times: user 0.05 s, sys: 0.00 s, total: 0.05 s
    Wall time: 0.05

GMP 表现稍好（预料之中，因为为 Sage 构建的 PARI 版本使用 GMP 进行整数运算）。

还可以使用 ``cputime`` 命令计时一组命令块，如下所示：

::

    sage: t = cputime()
    sage: a = int(1938)^int(99484)
    sage: b = 1938^99484
    sage: c = pari(1938)^pari(99484)
    sage: cputime(t)                       # somewhat random output
    0.64

.. skip

::

    sage: cputime?
    ...
        Return the time in CPU second since Sage started, or with optional
        argument t, return the time since time t.
        INPUT:
            t -- (optional) float, time in CPU seconds
        OUTPUT:
            float -- time in CPU seconds

``walltime`` 命令的行为与 ``cputime`` 命令类似，只是它计算的是挂钟时间。

我们也可以用 Sage 包含的计算机代数系统计算上面的幂。以下每种情况下，我们执行一个简单命令以启动该程序的服务器。
最相关的时间是挂钟时间。然而，如果挂钟时间和 CPU 时间之间存在显著差异，则可能表明存在值得优化的性能问题。

.. skip

::

    sage: time 1938^99484;
    CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
    Wall time: 0.01
    sage: gp(0)
    0
    sage: time g = gp('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: maxima(0)
    0
    sage: time g = maxima('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.30
    sage: kash(0)
    0
    sage: time g = kash('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: mathematica(0)
            0
    sage: time g = mathematica('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.03
    sage: maple(0)
    0
    sage: time g = maple('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.11
    sage: libgap(0)
    0
    sage: time g = libgap.eval('1938^99484;')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 1.02

注意，在这项测试中 GAP 和 Maxima 最慢（运行在 ``sage.math.washington.edu`` 机器上）。
由于 pexpect 接口的开销，将它们与最快的 Sage 相比可能不太公平。

其他 IPython 技巧
====================

如上文所述，Sage 使用 IPython 作为前端，因此你可以使用任何 IPython 的命令和功能。
你可以阅读
`完整的 IPython 文档 <http://ipython.scipy.org/moin/Documentation>`_ 。
下面是一些有趣的技巧 -- 在 IPython 中，这些被称为 "Magic 命令"：

- 如果你想输入一些复杂代码，可以使用 ``%edit`` （或 ``%ed`` 或 ``ed``）打开一个编辑器。
  在启动 Sage 之前，请确保 :envvar:`EDITOR` 环境变量设置为你喜欢的编辑器
  （通过在适当位置如 ``.profile`` 文件中放置 ``export EDITOR=/usr/bin/emacs`` 或
  ``export EDITOR=/usr/bin/vim`` 等）。在 Sage 提示符下执行 ``%edit`` 会打开指定的编辑器。
  然后在编辑器中你可以定义一个函数：

  .. CODE-BLOCK:: python

    def some_function(n):
        return n**2 + 3*n + 2

  保存并退出编辑器。在剩下的 Sage 会话期间，你可以使用 ``some_function``。
  如果你想修改它，可以在 Sage 提示符下输入 ``%edit some_function``。

- 如果你有一个计算，并且想修改其输出以便用于其他用途，可执行计算并输入 ``%rep``：
  这会将上一个命令的输出放置到 Sage 提示符，供你编辑。::

    sage: f(x) = cos(x)
    sage: f(x).derivative(x)
    -sin(x)

  此时如果你在 Sage 提示符下输入 ``%rep``, 你会得到一个新的 Sage 提示符，后面跟着 ``-sin(x)``, 光标在行尾。

要了解更多信息，请输入 ``%quickref`` 以获得 IPython 快速参考指南。
截止本文撰写时间（2011 年 4 月），Sage 使用的 IPython 版本为 0.9.1，
`Magic 命令文档 <http://ipython.org/ipython-doc/dev/interactive/tutorial.html#magic-functions>`_
可以在线访问。各种较为高级的 Magic 命令系统的内容记载在
`这里 <http://ipython.org/ipython-doc/stable/interactive/reference.html#magic-command-system>`_ 。


错误与异常
=====================

出现问题时，通常会看到 Python “异常”。Python 甚至会尝试给出引发异常的原因。
通常可以看到异常的名称，例如：:class:`NameError` 或 :class:`ValueError`
（详细异常列表请参见 Python 库参考 [PyLR]_ ）。例如：

.. skip

::

    sage: 3_2
    ------------------------------------------------------------
       File "<console>", line 1
         ZZ(3)_2
               ^
    SyntaxError: invalid ...

    sage: EllipticCurve([0,infinity])
    ------------------------------------------------------------
    Traceback (most recent call last):
    ...
    TypeError: Unable to coerce Infinity (<class 'sage...Infinity'>) to Rational

有时交互式调试器对理解问题很有用。可以使用 ``%pdb`` 切换它（默认是关闭的）。
如果打开调试器，出现异常时会出现提示符 ``ipdb>``。在调试器中，可以打印任意局部变量的状态，
并在执行栈中上下移动。例如：

.. skip

::

    sage: %pdb
    Automatic pdb calling has been turned ON
    sage: EllipticCurve([1,infinity])
    ---------------------------------------------------------------------------
    <class 'exceptions.TypeError'>             Traceback (most recent call last)
    ...

    ipdb>

在 ``ipdb>`` 提示符下输入 ``?`` 以获取调试器命令列表：

.. CODE-BLOCK:: text

    ipdb> ?

    Documented commands (type help <topic>):
    ========================================
    EOF    break  commands   debug    h       l     pdef   quit    tbreak
    a      bt     condition  disable  help    list  pdoc   r       u
    alias  c      cont       down     ignore  n     pinfo  return  unalias
    args   cl     continue   enable   j       next  pp     s       up
    b      clear  d          exit     jump    p     q      step    w
    whatis where

    Miscellaneous help topics:
    ==========================
    exec  pdb

    Undocumented commands:
    ======================
    retval  rv

输入 Ctrl-D 或 ``quit`` 返回 Sage。

.. _section-tabcompletion:

反向搜索与 Tab 补全
=================================

反向搜索：
输入命令的开头，然后按 ``Ctrl-p`` （或直接按上箭头键）查看以前输入的以该命令开头的命令行。
即使你完全退出 Sage 并稍后重新启动，这些功能仍然可以使用。也可以使用 ``Ctrl-r`` 通过历史记录进行反向搜索。
所有这些功能均使用 ``readline`` 软件包，可在大多数 Linux 版本中使用。

为了演示 Tab 补全，首先创建三维向量空间 :math:`V=\QQ^3` 如下：

::

    sage: V = VectorSpace(QQ,3)
    sage: V
    Vector space of dimension 3 over Rational Field

也可以使用如下更简洁的表示法：

::

    sage: V = QQ^3

然后可以很容易地使用 Tab 补全列出 :math:`V` 的所有成员函数。只需输入 ``V.``, 然后按键盘上的 :kbd:`Tab` 键:

.. skip

::

    sage: V.[tab key]
    V._VectorSpace_generic__base_field
    ...
    V.ambient_space
    V.base_field
    V.base_ring
    V.basis
    V.coordinates
    ...
    V.zero_vector

如果输入函数的前几个字母，然后按 :kbd:`Tab` 键，只会显示以这些字母开头的函数。

.. skip

::

    sage: V.i[tab key]
    V.is_ambient  V.is_dense    V.is_full     V.is_sparse

如果想知道某函数的作用，例如 coordinates 函数，
输入 ``V.coordinates?`` 来获取帮助或 ``V.coordinates??`` 查看源码，如下一节所述。


集成帮助系统
======================

Sage 拥有集成帮助系统。输入函数名后跟 ? 可以查看该函数的文档。

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates?
    Type:           instancemethod
    Base Class:     <class 'instancemethod'>
    String Form:    <bound method FreeModule_ambient_field.coordinates of Vector
    space of dimension 3 over Rational Field>
    Namespace:      Interactive
    File:           /home/was/s/local/lib/python2.4/site-packages/sage/modules/f
    ree_module.py
    Definition:     V.coordinates(self, v)
    Docstring:
        Write v in terms of the basis for self.

        Returns a list c such that if B is the basis for self, then

                sum c_i B_i = v.

        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: M = FreeModule(IntegerRing(), 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0-M1)
            [2, -1]

如上所示，输出告诉你对象的类型，定义它的文件，以及有用的函数描述及示例，
可以将这些示例粘贴到当前会话中。几乎所有这些示例都会定期自动测试，以确保它们正常工作并完全按照描述运行。

另一个非常符合 Sage 开源精神的功能是，如果 ``f`` 是一个 Python 函数，
那么输入 ``f??`` 会显示定义 ``f`` 的源代码。例如：

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates??
    Type:           instancemethod
    ...
    Source:
    def coordinates(self, v):
            """
            Write $v$ in terms of the basis for self.
            ...
            """
            return self.coordinate_vector(v).list()

这告诉我们 ``coordinates`` 函数所做的就是调用 ``coordinate_vector`` 函数并将结果转换为列表。
``coordinate_vector`` 函数做什么？

.. skip

::

    sage: V = QQ^3
    sage: V.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            ...
            return self.ambient_vector_space()(v)

``coordinate_vector`` 函数将其输入强制转化环绕空间，
其效果是以 :math:`V` 的形式计算 :math:`v` 的系数向量。
空间 :math:`V` 已经是环绕空间，因为它就是 :math:`\QQ^3`。
子空间也有 ``coordinate_vector`` 函数，它是不同的。我们创建一个子空间并看到：

.. skip

::

    sage: V = QQ^3; W = V.span_of_basis([V.0, V.1])
    sage: W.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            """
             ...
            """
            # First find the coordinates of v wrt echelon basis.
            w = self.echelon_coordinate_vector(v)
            # Next use transformation matrix from echelon basis to
            # user basis.
            T = self.echelon_to_user_matrix()
            return T.linear_combination_of_rows(w)

（如果你认为实现效率低下，请注册以帮助优化线性代数。）

你也可以输入 ``help(command_name)`` 或 ``help(class)`` 来获取给定类的帮助文档（类似 manpage ）。

.. skip

::

    sage: help(VectorSpace)
    Help on function VectorSpace in module sage.modules.free_module:

    VectorSpace(K, dimension_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
                with_basis='standard', dimension=None, basis_keys=None, **args)
    EXAMPLES:

    The base can be complicated, as long as it is a field.

    ::

        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x
         over Integer Ring
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
    --More--

当你输入 ``q`` 退出帮助系统时，你的会话内容将保持不变。
帮助列表不会使你的会话变得杂乱，而 ``function_name?`` 的输出有时会造成这种情况。
输入 ``help(module_name)`` 特别有用。例如，向量空间在 ``sage.modules.free_module`` 中定义，
输入 ``help(sage.modules.free_module)`` 即可获得有关整个模块的文档。
使用帮助查看文档时，可以通过输入 ``/`` 进行搜索，也可以通过输入 ``?`` 反向搜索。

保存和加载单个对象
=====================================

假设你计算出一个矩阵或更复杂的模符号空间，并希望将其保存以供日后使用。你要怎么办呢？
计算机代数系统采用多种方法来保存单个对象。


#. **保存游戏：** 仅支持保存和加载完整会话（如 GAP、Magma）。

#. **统一输入输出：** 使每个对象都以可读的方式打印（GP/PARI）。

#. **Eval:** 轻松在解释器中计算任意代码（如 Singular、PARI）。


由于 Sage 使用 Python，因此采用不同的方法，即每个对象都可以序列化，
转化为一个可以从中恢复该对象的字符串。这与 PARI 的统一输入输出方法精神相似，
只不过对象打印到屏幕的方式不会过于复杂。此外，保存和加载在大多数情况下是完全自动的，
不需要额外编程；这是 Python 的设计特性。

几乎所有 Sage 对象 x 都可以以压缩形式保存到磁盘，
使用 ``save(x, filename)`` （或在许多情况下 ``x.save(filename)``）。
要加载对象，使用 ``load(filename)``。

.. skip

::

    sage: A = MatrixSpace(QQ,3)(range(9))^2
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]
    sage: save(A, 'A')

现在你应该退出 Sage 并重新启动。然后便可以恢复 ``A``：

.. skip

::

    sage: A = load('A')
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]

可以使用同样的方法处理更复杂的对象，如椭圆曲线。缓存对象的所有数据都与对象一同保存。例如：

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: v = E.anlist(100000)              # takes a while
    sage: save(E, 'E')
    sage: quit

``E`` 的存储版占 153K 字节，因为它储存了前 100000 个 :math:`a_n`.

.. skip

::

    ~/tmp$ ls -l E.sobj
    -rw-r--r--  1 was was 153500 2006-01-28 19:23 E.sobj
    ~/tmp$ sage [...]
    sage: E = load('E')
    sage: v = E.anlist(100000)              # instant!

（在 Python 中，保存和加载使用 ``cPickle`` 模块实现。
具体来说，Sage 对象 ``x`` 可以通过 ``cPickle.dumps(x, 2)`` 保存。注意 ``2``！）

Sage 无法保存和加载某些其它计算机代数系统（例如 GAP、Singular、Maxima）创建的单个对象。
它们重新加载时状态显示为“无效 (invalid)”。在 GAP 中，虽然许多对象的打印方式可以重新构建，
但很多对象却不行，因此特意不允许从其打印表示进行重建。

.. skip

::

    sage: a = libgap(2)
    sage: a.save('a')
    sage: load('a')
    Traceback (most recent call last):
    ...
    ValueError: The session in which this object was defined is no longer
    running.

GP/PARI 对象可以保存和加载，因为它们的打印表示足以重构它们。

.. skip

::

    sage: a = gp(2)
    sage: a.save('a')
    sage: load('a')
    2

保存的对象稍后可以在不同架构或操作系统的计算机上重新加载，
例如，你可以在 32 位 OS X 上保存一个大矩阵，然后在 64 位 Linux 上重新加载它，
计算阶梯形式，然后再保存回去。此外，在许多情况下，即使在不同版本的 Sage 中也能加载对象，
只要该对象的代码没有太大差异。对象的所有属性，以及定义对象的类（但不包括源代码）都会被保存。
如果该类在新版本的 Sage 中不再存在，那么该对象就无法在新版本中重新加载。
但你可以在老版本中加载它，获取其对象字典（使用 ``x.__dict__``），保存该字典，并将其加载到新版本中。

保存为文本
--------------

你还可以将对象的 ASCII 文本表示保存到纯文本文件中，
只需以写入模式打开文件并写入对象的字符串表示即可（你也可以通过这种方式写入许多对象）。
写完对象后，关闭文件。

.. skip

::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: f = (x+y)^7
    sage: o = open('file.txt','w')
    sage: o.write(str(f))
    sage: o.close()

.. _section-save:

保存和加载完整会话
====================================

Sage 对于保存和加载完整会话有非常灵活的支持。

``save_session(sessionname)`` 命令将所有在当前会话中定义的变量保存为给定 ``sessionname`` 的字典。
（在少数情况下，如果某个变量不支持保存，则不会保存到字典。）生成的文件为 ``.sobj`` 文件，
可以像保存的其它对象一样加载。加载会话保存的对象时，会得到一个字典，字典的键为变量名，值为对象。

可以使用 ``load_session(sessionname)`` 命令将 ``sessionname`` 中定义的变量加载到当前会话。
注意，这不会清除当前会话中已经定义的变量；而是合并两个会话。

首先启动 Sage 并定义一些变量。

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: M = ModularSymbols(37)
    sage: a = 389
    sage: t = M.T(2003).matrix(); t.charpoly().factor()
     _4 = (x - 2004) * (x - 12)^2 * (x + 54)^2

接下来保存会话，将上面定义的每个变量保存至文件。然后查看文件，大小约为 3K。

.. skip

::

    sage: save_session('misc')
    Saving a
    Saving M
    Saving t
    Saving E
    sage: quit
    was@form:~/tmp$ ls -l misc.sobj
    -rw-r--r--  1 was was 2979 2006-01-28 19:47 misc.sobj

最后重新启动 Sage，定义一个额外的变量，并加载保存的会话。

.. skip

::

    sage: b = 19
    sage: load_session('misc')
    Loading a
    Loading M
    Loading E
    Loading t

每个保存的变量再次可用。此外，变量 ``b`` 没有被覆盖。

.. skip

::

    sage: M
    Full Modular Symbols space for Gamma_0(37) of weight 2 with sign 0
    and dimension 5 over Rational Field
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational
    Field
    sage: b
    19
    sage: a
    389

