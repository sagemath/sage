***********************
Sage, LaTeX 及其朋友们
***********************

Sage 与 TeX 的 LaTeX 方言之间存在着密切的协同关系。
本节旨在介绍各种交互方式，从最基本的开始，然后介绍一些不常见的用法。

基本使用
=========

Sage 中的每个“对象”都必须有 LaTeX 表示。你可以通过执行 ``latex(foo)`` 来获取这种表示，
其中 ``foo`` 是 Sage 中的某个对象。输出是一个字符串，当在 TeX 的数学模式中使用时
（例如，包围在一对单美元符号之间），该字符串应该能够准确地呈现 ``foo``。以下是一些示例。 ::

    sage: var('z')
    z
    sage: latex(z^12)
    z^{12}
    sage: latex(sqrt(z^2 + 1/2))
    \sqrt{z^{2} + \frac{1}{2}}
    sage: latex('a string')
    \text{\texttt{a{ }string}}
    sage: latex(QQ)
    \Bold{Q}
    sage: latex(ZZ['x'])
    \Bold{Z}[x]
    sage: latex(matrix(QQ, 2, 3, [[2,4,6],[-1,-1,-1]]))
    \left(\begin{array}{rrr}
    2 & 4 & 6 \\
    -1 & -1 & -1
    \end{array}\right)

通过这种方式，Sage 可以有效地用于构建 LaTeX 文档的各个部分：
在 Sage 中创建或计算一个对象 ``foo``，对该对象执行 ``latex(foo)``，
然后将 LaTeX 字符串剪切并粘贴到你的文档中。

命令 ``view(foo)`` 会显示对象 ``foo`` 的渲染后的 LaTeX 表示。
在后台，该命令会运行 ``latex(foo)`` 并将 LaTeX 字符串合并到一个简单的 LaTeX 文档中，
用系统范围内的 TeX 安装处理该文档，然后调用合适的查看器来显示输出。

在 Jupyter Notebook 中，你可以自动看到输入命令输出的渲染 LaTeX 表示。
你可以通过执行 ``%display latex`` 来启动自动渲染（并通过执行 ``%display plain`` 停止）。

.. ONLY:: html and feature_jupyter_sphinx

    因此，在 Jupyter notebook 中，你得到

    .. JUPYTER-EXECUTE::

        %display latex
        var('z')
        z^12

    .. JUPYTER-EXECUTE::

        sqrt(z^2 + 1/2)

    .. JUPYTER-EXECUTE::

        'a string'

    .. JUPYTER-EXECUTE::

        QQ

    .. JUPYTER-EXECUTE::

        ZZ['x']

    .. JUPYTER-EXECUTE::

        matrix(QQ, 2, 3, [[2,4,6],[-1,-1,-1]])

    .. JUPYTER-EXECUTE::

        %display plain

Jupyter Notebook 使用 `MathJax <http://www.mathjax.org>`_ 在网页浏览器中清晰地渲染数学内容。
MathJax 是一个开源的 JavaScript 数学显示引擎，可以在所有现代浏览器中使用。
它能够渲染大部分 LaTex，但并不支持完整的 LaTeX，是 LaTex 的子集。
它不支持复杂表格、分段或文档管理，因为它主要用于准确渲染 LaTeX 数学片段。

在 Jupyter Notebook 中自动 LaTeX 渲染（启用 ``%display latex``）
是通过 :class:`sage.misc.html.MathJax` 类内部实现的。
该类的对象将 Sage 对象通过 ``latex()`` 转换为 MathJax 需要的 HTML 形式，然后将其包装在 HTML 中。  ::

    sage: from sage.misc.html import MathJax
    sage: mj = MathJax()
    sage: var('z')
    z
    sage: mj(z^12)
    <html>\[z^{12}\]</html>
    sage: mj(sqrt(z^2 + 1/2))
    <html>\[\sqrt{z^{2} + \frac{1}{2}}\]</html>
    sage: mj('a string')
    <html>\[\verb|a|\verb| |\verb|string|\]</html>
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}\]</html>
    sage: mj(ZZ['x'])
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}[x]\]</html>
    sage: mj(matrix(QQ, 2, 3, [[2,4,6],[-1,-1,-1]]))
    <html>\[\left(\begin{array}{rrr}
    2 & 4 & 6 \\
    -1 & -1 & -1
    \end{array}\right)\]</html>

如果你需要了解 Sage 对象的 LaTeX 渲染，那么了解这一点很有用。


.. _sec-custom-generation:

自定义 LaTeX 生成
============================

有几种方法可以自定义由 ``latex()`` 命令生成的实际 LaTeX 代码。
预定义对象 ``latex`` 包含多个方法，可以通过输入 ``latex.`` （注意这里有一个点）后按 :kbd:`Tab` 键来列出这些方法。

``latex.matrix_delimiters`` 方法是一个很好的例子。
它可以用来更改矩阵周围的符号 -- 大括号、方括号、花括号、竖线。
不强制执行任何样式，你可以随意混合搭配。
注意，LaTeX 所需的反斜杠在 Python 字符串中需要额外加一个斜杠以便正确转义。  ::

    sage: A = matrix(ZZ, 2, 2, range(4))
    sage: latex(A)
    \left(\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right)
    sage: latex.matrix_delimiters(left='[', right=']')
    sage: latex(A)
    \left[\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right]
    sage: latex.matrix_delimiters(left='\\{', right='\\}')
    sage: latex(A)
    \left\{\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right\}

``latex.vector_delimiters`` 方法的工作原理与之类似。

常见环和域（整数、有理数、实数等）的排版方式可以通过 ``latex.blackboard_bold`` 方法来控制。
这些集合默认以粗体排版，但有时可以选择以双重划线格式书写，如某些书面作品所做的那样。
这可以通过重新定义 Sage 内置的 ``\Bold{}`` 宏来实现。 ::

    sage: latex(QQ)
    \Bold{Q}
    sage: from sage.misc.html import MathJax
    sage: mj = MathJax()
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}\]</html>
    sage: latex.blackboard_bold(True)
    sage: mj(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbb{#1}}\Bold{Q}\]</html>
    sage: latex.blackboard_bold(False)

.. ONLY:: html

    在 Jupyter notebook 中，

    .. JUPYTER-EXECUTE::

        %display latex
        QQ

    .. JUPYTER-EXECUTE::

        latex.blackboard_bold(True)
        QQ

    .. JUPYTER-EXECUTE::

        latex.blackboard_bold(False)
        %display plain

可以通过加入新的宏来利用 LaTeX 的可扩展性。可以添加单个宏，以便在 MathJax 解释 LaTeX 片段时使用。 ::

    sage: latex.add_macro(r"\newcommand{\sqrt}[1]{(#1)^\frac{1}{2}}")
    sage: latex.extra_macros()
    '\\newcommand{\\sqrt}[1]{(#1)^\\frac{1}{2}}'
    sage: var('x y')
    (x, y)
    sage: latex(sqrt(x+y))
    \sqrt{x + y}
    sage: from sage.misc.html import MathJax
    sage: mj = MathJax()
    sage: mj(sqrt(x + y))
    <html>\[\newcommand{\sqrt}[1]{(#1)^\frac{1}{2}}\sqrt{x + y}\]</html>
    sage: latex.extra_macros('')

.. ONLY:: html

    在 Jupyter notebook 中，

    .. JUPYTER-EXECUTE::

        %display latex
        var('x y')
        sqrt(x + y)

    .. JUPYTER-EXECUTE::

        latex.add_macro(r"\newcommand{\sqrt}[1]{(#1)^\frac{1}{2}}")
        sqrt(x + y)

    .. JUPYTER-EXECUTE::

        latex.extra_macros('')
        %display plain


.. _sec-custom-processing:

自定义 LaTeX 处理
============================

系统范围内的 TeX 被调用来处理完整的 LaTeX 文档，例如，当你 ``view(foo)`` 时，
其中 ``foo`` 是一个复杂的 Sage 对象，太复杂以至于 ``MathJax`` 无法处理。
命令 ``latex_extra_preamble`` 用于构建完整 LaTeX 文档的导言部分，下面将展示如何完成这项工作。
如往常一样，请注意 Python 字符串中需要双反斜杠。 ::

    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: from sage.misc.latex import latex_extra_preamble
    sage: print(latex_extra_preamble())
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: print(latex_extra_preamble())
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    \newcommand{\foo}{bar}

同样，对于更大或更复杂的 LaTeX 表达式，可以将包（或其他任意内容）添加到 LaTeX 文件的导言部分。
任意内容都可以通过 ``latex.add_to_preamble`` 命令加入导言部分，
专用命令 ``latex.add_package_to_preamble_if_available`` 会首先检查某个包是否实际存在，
然后尝试将其添加到导言部分。

这里我们将几何包添加到导言部分并用它来设置 TeX 将在页面上使用的区域尺寸（有效地设置边距）。
如往常一样，请注意 Python 字符串中需要双反斜杠。 ::

    sage: from sage.misc.latex import latex_extra_preamble
    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: latex.add_to_preamble('\\usepackage{geometry}')
    sage: latex.add_to_preamble('\\geometry{letterpaper,total={8in,10in}}')
    sage: latex.extra_preamble()
    '\\usepackage{geometry}\\geometry{letterpaper,total={8in,10in}}'
    sage: print(latex_extra_preamble())
    \usepackage{geometry}\geometry{letterpaper,total={8in,10in}}
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}

可以通过检查其存在性来添加特定包，以下示例展示了这种情况。作为示例，我们将尝试向导言部分添加一个可能不存在的包。 ::

    sage: latex.extra_preamble('')
    sage: latex.extra_preamble()
    ''
    sage: latex.add_to_preamble('\\usepackage{foo-bar-unchecked}')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'
    sage: latex.add_package_to_preamble_if_available('foo-bar-checked')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'

使用哪种 TeX 方言，以及输出和相关查看器的性质，也可以定制。

.. NOTE::

    Sage 几乎包括了构建和使用 Sage 所需的一切，但一个重要的例外是 TeX 本身。
    因此，在以下情况下，你需要安装完整的 TeX 系统以及一些相关的转换工具。
    许多版本的 Linux 都有基于 TeXLive 的软件包，macOS 有 MacTeX，Windows 有 MiKTeX。

可以使用 ``latex.engine()`` 命令控制是否使用系统范围内的 ``latex``, ``pdflatex`` 或 ``xelatex`` 可执行文件。
当调用 ``view()`` 并且引擎设置为 ``latex`` 时，会生成一个 dvi 文件，Sage 会使用 dvi 查看器（如 xdvi）来显示结果。
相比之下，当引擎设置为 ``pdflatex`` 时，调用 ``view()`` 会生成 PDF 文件，
并且 Sage 会调用系统的 PDF 文件查看工具（如 acrobat, okular, evince 等）。

对于使用这些工具的练习，有一些预先打包好的示例。
要使用这些示例，需要导入 ``sage.misc.latex.latex_examples`` 对象，
这是 :class:`sage.misc.latex.LatexExamples` 类的一个实例，如下所示。
目前该类有交换图、组合图、扭结理论和 pstricks 的示例，分别使用以下包：xy，tkz-graph，xypic，pstricks。
导入后，对 ``latex_examples`` 使用 tab 补全查看内置示例。
调用每个示例会返回一些关于如何正确呈现该示例的说明。要实际查看示例，需要使用 ``view(foo)`` （导言部分、引擎等均设置正确）。 ::

    sage: from sage.misc.latex import latex_examples
    sage: foo = latex_examples.diagram()
    sage: foo
    LaTeX example for testing display of a commutative diagram produced
    by xypic.
    <BLANKLINE>
    To use, try to view this object -- it will not work.  Now try
    'latex.add_to_preamble("\\usepackage[matrix,arrow,curve,cmtip]{xy}")',
    and try viewing again. You should get a picture (a part of the diagram arising
    from a filtered chain complex).

为了展示如何处理复杂的 LaTeX 表达式，让我们看一下使用 ``tkz-graph`` LaTeX 包的组合图示例。

.. NOTE::

    ``tkz-graph`` LaTeX 包建立在 ``pgf`` 库的 ``tikz`` 前端之上。
    渲染组合图需要 ``pgf`` 库以及文件 ``tkz-graph.sty`` 和 ``tkz-berge.sty``。
    它们很可能已经是系统范围内 TeX 安装的一部分。即使不是，也应当很容易找到安装指南。

首先，我们通过将相关包添加到 LaTeX 文档的导言部分来确保它们被包含在内。 ::

    sage: latex.extra_preamble('\\usepackage{tikz}\n\\usepackage{tkz-graph}\n'
    ....:                      '\\usepackage{tkz-berge}\n\\usetikzlibrary{arrows,shapes}')

当使用 dvi 文件作为中间格式时，图形无法正确生成，因此最好将 LaTeX 引擎设置为 ``pdflatex`` 可执行文件。 ::

    sage: latex.engine('pdflatex')

此时，像 ``view(graphs.CompleteGraph(4))`` 这样的命令应该生成一个带有完整图 `K_4` 适当图像的 PDF。

实际上，可以省略前面的步骤，因为导言部分会自动正确设置，并且 ``pdflatex`` 是 Sage 的默认 LaTeX 引擎。
重新启动 Sage 后再次尝试该命令。

注意，通过 ``tkz-graph`` 有多种选项可以影响 LaTeX 中图形的呈现方式，这超出了本节的范围。
请参阅参考手册 :ref:`sage.graphs.graph_latex` 章节获取指令和详细信息。


SageTeX
=======

SageTeX 是一个可以进一步集成 TeX 和 Sage 的程序。
它是一组 TeX 宏，允许 LaTeX 文档包含指令，让 Sage 计算各种对象并使用 ``latex()`` 格式化对象。
更多信息请参见 :ref:`sec-sagetex`。
