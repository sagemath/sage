.. _sec-sagetex:

*************
使用 SageTeX
*************

SageTeX 包允许你将 Sage 计算结果嵌入到 LaTeX 文档中。
要使用它，需要先“安装”它（请参阅 :ref:`sec-sagetex_install`）。

示例
----------

以下是一个非常简短的 SageTeX 使用示例。
完整文档可以在 :file:`SAGE_ROOT/venv/share/doc/sagetex` 中找到，
其中 ``SAGE_ROOT`` 是 Sage 安装目录。该目录包含文档和示例文件。
请参阅 :file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex` 以获取一些可能有用的 Python 脚本。

想要了解 SageTeX 的工作原理，请按照 SageTeX 的安装说明（在 :ref:`sec-sagetex_install` 中）操作，
并将以下文本复制到一个名为 ``st_example.tex`` 的文件中：

.. warning::

  如果你在“实时”帮助中查看此内容，下面的文本会有几个未知控制序列的错误。
  请使用静态版查看正确的文本。

.. code-block:: latex

    \documentclass{article}
    \usepackage{sagetex}

    \begin{document}

    Using Sage\TeX, one can use Sage to compute things and put them into
    your \LaTeX{} document. For example, there are
    $\sage{number_of_partitions(1269)}$ integer partitions of $1269$.
    You don't need to compute the number yourself, or even cut and paste
    it from somewhere.

    Here's some Sage code:

    \begin{sageblock}
        f(x) = exp(x) * sin(2*x)
    \end{sageblock}

    The second derivative of $f$ is

    \[
      \frac{\mathrm{d}^{2}}{\mathrm{d}x^{2}} \sage{f(x)} =
      \sage{diff(f, x, 2)(x)}.
    \]

    Here's a plot of $f$ from $-1$ to $1$:

    \sageplot{plot(f, -1, 1)}

    \end{document}

像往常一样在 ``st_example.tex`` 上运行 LaTeX。请注意 LaTeX 会有一些警告，其中包括：

.. CODE-BLOCK:: text

    Package sagetex Warning: Graphics file
    sage-plots-for-st_example.tex/plot-0.eps on page 1 does not exist. Plot
    command is on input line 25.

    Package sagetex Warning: There were undefined Sage formulas and/or
    plots. Run Sage on st_example.sagetex.sage, and then run LaTeX on
    st_example.tex again.

请注意，除了 LaTeX 产生的常规文件集合外，还有一个名为 ``st_example.sagetex.sage`` 的文件。
这是在 ``st_example.tex`` 上运行 LaTeX 时生成的 Sage 脚本。
警告信息告诉你在 ``st_example.sagetex.sage`` 上运行 Sage，请听从建议并进行操作。
它会告诉你再次在 ``st_example.tex`` 上运行 LaTeX，但在此之前，
请注意新文件 ``st_example.sagetex.sout`` 已被创建。该文件包含 Sage 计算结果，
可供 LaTeX 插入到你的文本中。还创建了一个包含 EPS 文件的新目录。再次运行 LaTeX，
你会看到 Sage 计算和绘图的所有内容已包含在你的文档中。

上面使用的各种宏应该很容易理解。``sageblock`` 环境按原样排版你的代码，
并在运行 Sage 时执行代码。当你执行 ``\sage{foo}`` 时，
插入到文档中的结果就是在 Sage 内部运行 ``latex(foo)`` 得到的结果。
绘图命令稍微复杂一些，但在最简单形式下，``\sageplot{foo}`` 插入的是由 ``foo.save('filename.eps')`` 得到的图像。

一般来说，操作步骤是：

    - 在 .tex 文件上运行 LaTeX；
    - 在生成的 .sage 文件上运行 Sage；
    - 再次运行 LaTeX。

如果文档中没有更改任何 Sage 命令，则可以省略运行 Sage。

SageTeX 还有很多内容，由于 Sage 和 LaTeX 都是复杂且强大的工具，
建议阅读 SageTeX 的文档 :file:`SAGE_ROOT/venv/share/doc/sagetex`。

.. _sec-sagetex_install:

让 TeX 识别 SageTeX
-------------------------

Sage 基本上是自包含的，但某些部分需要进行一些干预才能正常工作。SageTeX 就是其中之一。

SageTeX 包允许在 LaTeX 文档中嵌入来自 Sage 的计算和绘图。
Sage 中默认安装了 SageTeX，但要在 LaTeX 文档中使用 SageTeX，你需要先让 TeX 识别它。

关键在于 TeX 需要能够找到 sagetex.sty，
该文件位于 :file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/`，
其中 ``SAGE_ROOT`` 是你构建或安装 Sage 的目录。如果 TeX 能找到 ``sagetex.sty``，
那么 SageTeX 就可以工作。有几种方法可以实现这一点。

- 第一种方法，也是最简单的方法是将 ``sagetex.sty`` 复制到与 LaTeX 文档相同的目录中。
  在排版文档时，总会搜索当前目录，因此这种方法始终有效。

  但这种方法有两个小问题：首先，会在计算机上产生很多不必要的 ``sagetex.sty`` 拷贝。
  其次，更严重的问题是，如果升级 Sage 并获得新版本的 SageTeX，
  Python 代码和 SageTeX 的 LaTeX 代码可能不再匹配，从而导致错误。

- 第二种方法是使用 ``TEXMFLOCAL`` 环境变量。如果你使用的是 bash shell，可以这样做：

  .. CODE-BLOCK:: shell-session

      $ export TEXMFLOCAL=SAGE_ROOT/venv/share/texmf
      $ mktexlsr       # update kpathsea ls-R databases

  其中 ``SAGE_ROOT`` 是 Sage 安装位置。
  之后，TeX 和相关程序将找到 SageTeX 样式文件。如果你想使这个更改持续生效，
  可以将上述第一行添加到 ``.bashrc`` 文件中。如果你使用的是不同的 shell，
  可能需要调整以上命令从而让环境变量可被识别；请查阅所用 shell 的文档以了解如何操作。

  如果你移动了 Sage 的安装目录或在新目录中安装了新版本，
  需要用新的 ``SAGE_ROOT`` 更新上述命令。

- 让 TeX 识别 ``sagetex.sty`` 的第三种（也是最佳的）方法，
  是将该文件复制到主目录中的一个方便的位置。
  大多数 TeX 发行版会自动搜索主目录中的 ``texmf`` 目录以寻找包。
  要确切了解这个目录的位置，请在命令行种执行以下操作：

  .. CODE-BLOCK:: shell-session

      $ kpsewhich -var-value=TEXMFHOME

  这将打印出一个目录，例如 ``/home/drake/texmf`` 或 ``/Users/drake/Library/texmf``。
  使用如下命令将 :file:`SAGE_ROOT/venv/share/texmf/` 中的 ``tex/`` 目录复制到主目录的 ``texmf`` 目录：

  .. CODE-BLOCK:: shell-session

      $ cp -R SAGE_ROOT/venv/share/texmf/tex TEXMFHOME

  其中 ``SAGE_ROOT`` 仍然是 Sage 的安装位置，``TEXMFHOME`` 是 ``kpsewhich`` 命令的结果。

  如果你升级了 Sage 并发现 SageTeX 无法工作，
  可以简单地重复上述步骤以确保 SageTeX 的 Sage 部分和 TeX 部分再次同步。

.. _sagetex_installation_multiuser:

- 对于多用户系统上的安装，只需适当修改上述指令，将 ``sagetex.sty`` 复制到系统范围的 TeX 目录中。
  最好的选择可能是使用以下结果，而不是 ``TEXMFHOME`` 目录：

  .. CODE-BLOCK:: shell-session

      $ kpsewhich -var-value=TEXMFLOCAL

  这很可能会产生类似于 ``/usr/local/share/texmf`` 的结果。
  按照上述方式将 ``tex`` 目录复制到 ``TEXMFLOCAL`` 目录中。
  现在需要通过运行以下命令更新 TeX 的包数据库：

  .. CODE-BLOCK:: shell-session

      $ texhash TEXMFLOCAL

  以 root 身份，适当替换 ``TEXMFLOCAL``。
  现在系统中所有用户都可以访问 LaTeX 包，如果他们也能运行 Sage，他们就可以使用 SageTeX。

.. warning::

  确保 LaTeX 在排版文档时使用的 ``sagetex.sty`` 文件与 SageTeX 使用的版本匹配，
  这一点至关重要。如果你升级了 Sage，应该删除所有旧版本的 ``sagetex.sty``。

  由于此问题，我们建议将 SageTeX 文件复制到主目录的 texmf 目录中（上述第 3 种方法）。
  这样，升级 Sage 时，仅需做一件事（复制目录）即可确保 SageTeX 正常工作。

SageTeX 文档
---------------------

虽然这不严格属于安装的一部分，但值得在此提及的是，
SageTeX 的文档维护在 :file:`SAGE_ROOT/venv/share/doc/sagetex/sagetex.pdf`。
同一目录中还有一个示例文件 -- 请参见 ``example.tex`` 和 ``example.pdf``，
这是使用 LaTeX 和 Sage 对该文件进行排版的预生成结果。
你也可以从 `SageTeX 页面 <https://github.com/sagemath/sagetex>`_ 获取这些文件。

SageTeX 与 TeXLive
-------------------

一个潜在的令人困惑的问题是流行的 TeX 发行版 `TeXLive <http://www.tug.org/texlive/>`_ 包含 SageTeX。
虽然看起来很方便，但对于 SageTeX 而言，确保 Sage 部分和 LaTeX 部分同步是非常重要的 -- 在这种情况下，
这就成为了一个问题，因为由操作系统发行版或软件包管理器提供的 TeXLive 可能与官方 TeXLive 分发版本不同步，
而后者也可能与当前的 SageTeX 版本不同步。

因此，*强烈建议* 你始终按照上面的说明，从 Sage 安装 SageTeX 的 LaTeX 部分。
上述说明将确保 SageTeX 的两个部分兼容并正常工作。
