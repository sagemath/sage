************
介绍
************

完成本教程最多需要 3-4 小时。你可以阅读本教程的 HTML 或 PDF 版本，
或者在 Sage Notebook 中点击 ``Help``，然后点击 ``Tutorial`` 以交互方式在 Sage 中完成教程。

虽然 Sage 的大部分是使用 Python 实现的，但阅读本教程并不需要 Python 背景。
可能你会在某个时点希望学习 Python（一门非常有趣的语言！），有很多优秀的免费资源可以帮助你：
Python 初学者指南 [PyB]_ 列出了许多选择。如果你只是想快速试用 Sage，那么本教程是很好的起点。例如：

::

    sage: 2 + 2
    4
    sage: factor(-2007)
    -1 * 3^2 * 223

    sage: A = matrix(4,4, range(16)); A
    [ 0  1  2  3]
    [ 4  5  6  7]
    [ 8  9 10 11]
    [12 13 14 15]

    sage: factor(A.charpoly())
    x^2 * (x^2 - 30*x - 80)

    sage: m = matrix(ZZ,2, range(4))
    sage: m[0,0] = m[0,0] - 3
    sage: m
    [-3  1]
    [ 2  3]

    sage: E = EllipticCurve([1,2,3,4,5]);
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
    over Rational Field
    sage: E.anlist(10)
    [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    sage: E.rank()
    1

    sage: k = 1/(sqrt(3)*I + 3/4 + sqrt(73)*5/9); k
    36/(20*sqrt(73) + 36*I*sqrt(3) + 27)
    sage: N(k)
    0.165495678130644 - 0.0521492082074256*I
    sage: N(k,30)      # 30 "bits"
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{36}{20 \, \sqrt{73} + 36 i \, \sqrt{3} + 27}

.. _installation:

安装
============

如果你的电脑没有安装 Sage，只是想尝试一些命令，可以在 http://sagecell.sagemath.org 上在线使用。

请参阅 Sage 主页 [SA]_ 文档中的安装指南，了解如何在你的电脑上安装 Sage。以下是一些简要说明。

#. Sage 下载文件附带所有所需组件。换句话说，虽然 Sage 使用 Python、IPython、PARI、GAP、Singular、Maxima、NTL、GMP 等，
   你不需要单独安装它们，因为它们已经包含在 Sage 发行版中。
   但是，要使用某些 Sage 功能，例如 Macaulay 或 KASH，你必须确保电脑已经安装了相关程序。

#. Sage 的预编译二进制版本（可以在 Sage 官网上找到）可能比源代码版本更容易和更快安装。只需解压文件并运行 ``sage``。

#. 如果你想使用 SageTeX 包（允许你将 Sage 计算结果嵌入到 LaTeX 文件中），
   你需要让 TeX 发行版识别 SageTeX。请参阅 `Sage 安装指南 <http://doc.sagemath.org/html/en/>`_
   中的“让 TeX 识别 SageTeX”章节（这个链接 `<../installation/index.html>`_ 会为你打开安装指南）。
   其实非常简单；你只需要设置一个环境变量或复制一个文件到 TeX 的搜索目录中。

   如何使用 SageTeX 的文档位于 ``$SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/``，
   其中 "``$SAGE_ROOT``" 指的是 Sage 的安装目录，例如 ``/opt/sage-9.6``。


使用 Sage 的方法
================

Sage 可以通过多种方式使用：


-  **Notebook 图形界面：** 运行 ``sage -n jupyter``; 请参阅
   `Jupyter 在线文档 <https://jupyter-notebook.readthedocs.io/en/latest/notebook.html>`_,

-  **交互式 Shell：** 请参阅 :ref:`chapter-interactive_shell`,

-  **编写程序：** 在 Sage 中编写解释和编译的程序（请参阅 :ref:`section-loadattach` 和 :ref:`section-compile`)

-  **编写脚本：** 编写使用 Sage 库的独立 Python 脚本（请参阅 :ref:`section-standalone`).


Sage 的长期目标
=======================

-  **实用：** Sage 的目标受众包括学习数学的学生（从高中到研究生）、教师和研究数学家。
   旨在提供可以用于探索和实验代数、几何、数论、微积分、数值计算等数学构造的软件。
   Sage 能够帮助用户更方便地进行数学对象的交互实验。

-  **高效：** Sage 追求快速。它使用高度优化的成熟软件，如 GMP、PARI、GAP 和 NTL，因此在某些操作上非常快速。

-  **免费开源：** 源代码必须免费提供且可读，用户可以了解系统的实际运作，并能够更容易地进行扩展。
   正如数学家通过仔细阅读或浏览证明来深入理解定理一样，用户应该能够通过阅读带有文档的源代码来理解计算过程。
   如果在发表的论文中使用 Sage 进行计算，读者将始终可以免费访问 Sage 及其所有源代码，
   你甚至可以归档和重新分发你使用的 Sage 版本。

-  **易于编译：** Sage 应该易于从源代码编译，适用于 Linux、OS X 和 Windows 用户，这使得用户修改系统更加灵活。

-  **协作：** 提供与大多数其他计算机代数系统的强大接口，
   包括 PARI、GAP、Singular、Maxima、KASH、Magma、Maple 和 Mathematica。
   Sage 旨在统一和扩展现有数学软件。

-  **文档齐全：** 提供教程、编程指南、参考手册和操作指南，包含大量示例和背景数学讨论。

-  **可扩展：** 能够定义新的数据类型或从内置类型派生，并能够使用多种编程语言编写的代码。

-  **用户友好：** 功能易于理解，文档和源代码易于查看，并且提供高水平的用户支持。

