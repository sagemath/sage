***********
编程
***********

.. _section-loadattach:

加载和附加 Sage 文件
================================

接下来我们说明如何将写在单独文件中的程序加载到 Sage 中。
创建一个名为 ``example.sage`` 的文件，并写入以下内容：

.. CODE-BLOCK:: python

    print("Hello World")
    print(2^3)

你可以使用 ``load`` 命令读取并执行 ``example.sage`` 文件。

.. skip

::

    sage: load("example.sage")
    Hello World
    8

你也可以使用 ``attach`` 命令将 Sage 文件附加到运行的会话中：

.. skip

::

    sage: attach("example.sage")
    Hello World
    8

现在如果你修改 ``example.sage`` 文件并在 Sage 中输入一个空行（即按下回车键），
那么 ``example.sage`` 的内容将会自动重新加载到 Sage 中。

特别是，``attach`` 命令会在文件更改时自动重新加载文件，这在调试代码时非常方便，
而 ``load`` 命令仅加载文件一次。

当 Sage 加载 ``example.sage`` 时，它会将其转换为 Python，然后由 Python 解释器执行。
此转换非常简单；它主要是将整型字面量包装在 ``Integer()`` 中，
将浮点型字面量包装在 ``RealNumber()`` 中，
将 ``^`` 替换为 ``**``，并将例如 ``R.2`` 替换为 ``R.gen(2)``。
转换后的 ``example.sage`` 版本包含在与 ``example.sage`` 相同的目录中，
名为 ``example.sage.py``。该文件包含以下代码：

.. CODE-BLOCK:: python

    print("Hello World")
    print(Integer(2)**Integer(3))

整型字面量被包装，``^`` 被替换为 ``**``。
（在 Python 中，``^`` 表示“异或”，而 ``**`` 表示“幂运算”。）

（这种预解析由 ``sage/misc/interpreter.py`` 模块实现。）

只要有换行符来创建新块（在文件中则无需如此），你就可以将多行缩进代码粘贴到 Sage 中。
然而，更好的方式是将这些代码保存到文件中，并如上所述使用 ``attach`` 命令来加载。


.. _section-compile:

创建编译代码
======================

速度在数学计算中至关重要，因为更快的计算可以大大提高效率。
尽管 Python 是一种非常方便的高级语言，但如果使用静态类型的编译型语言实现某些计算，
其速度可以比用 Python 实现快几个数量级。如果 Sage 完全用 Python 编写，
那么在某些方面速度会过于缓慢。为了应对这种情况，Sage 支持一种编译“版本”的 Python，称为 Cython ([Cyt]_ 和 [Pyr]_)。
Cython 类似于 Python 和 C 语言。大多数 Python 结构，包括列表推导式、条件表达式、
类似 ``+=`` 这样的代码都支持；你还可以导入其他 Python 模块中编写的代码。
此外，你还可以声明任意的 C 变量，并直接调用任意的 C 库函数。生成的代码会转换为 C，并使用 C 编译器进行编译。

为了创建你自己的编译 Sage 代码，请将文件命名为 ``.spyx`` 扩展名（而非 ``.sage``）。
如果使用命令行界面，你可以像处理解释代码一样附加和加载编译代码（目前，Notebook 界面不支持附加和加载 Cython 代码）。
实际编译是在“后台”完成的，你无需进行任何显式操作。编译后的共享对象库存储在 ``$HOME/.sage/temp/hostname/pid/spyx`` 中。
这些文件将在退出 Sage 时删除。

Sage 预解析不适用于 spyx 文件，例如，``1/3`` 在 spyx 文件中结果为 0，
而不是有理数 :math:`1/3`。如果 ``foo`` 是 Sage 库中的一个函数，要想在 spyx 文件中使用它，
请导入 ``sage.all`` 并使用 ``sage.all.foo``。

.. CODE-BLOCK:: python

    import sage.all
    def foo(n):
        return sage.all.factorial(n)

访问单独文件中的 C 函数
---------------------------------------

访问定义在单独 \*.c 文件中的 C 函数也很容易。
以下是一个示例。在同一目录下创建文件 ``test.c`` 和 ``test.spyx``，内容如下：

纯 C 代码：``test.c``

.. CODE-BLOCK:: c

    int add_one(int n) {
      return n + 1;
    }

Cython 代码：``test.spyx``:

.. CODE-BLOCK:: cython

    cdef extern from "test.c":
        int add_one(int n)

    def test(n):
        return add_one(n)

然后进行以下操作：

.. skip

::

    sage: attach("test.spyx")
    Compiling (...)/test.spyx...
    sage: test(10)
    11

如果需要额外的库 ``foo`` 来编译从 Cython 文件生成的 C 代码，
在 Cython 源代码中添加 ``clib foo``。
类似地，可以使用声明 ``cfile bar`` 将额外的 C 文件 ``bar`` 包含在编译中。

.. _section-standalone:

独立 Python/Sage 脚本
==============================

以下独立 Sage 脚本可以分解整数、多项式等：

.. CODE-BLOCK:: python

    #!/usr/bin/env sage

    import sys

    if len(sys.argv) != 2:
        print("Usage: %s <n>" % sys.argv[0])
        print("Outputs the prime factorization of n.")
        sys.exit(1)

    print(factor(sage_eval(sys.argv[1])))

为了使用此脚本，``SAGE_ROOT`` 必须包含在 PATH 中。如果将上述脚本命名为 ``factor``，则以下是使用示例：

.. CODE-BLOCK:: shell-session

    $ ./factor 2006
    2 * 17 * 59

数据类型
==========

在 Sage 中，每个对象都有一个明确的类型。Python 有各种基本内置类型，
而 Sage 库还增加了更多类型。Python 内置类型包括字符串、列表、元组、整型和浮点型等，如下所示：

::

    sage: s = "sage"; type(s)
    <... 'str'>
    sage: s = 'sage'; type(s)      # you can use either single or double quotes
    <... 'str'>
    sage: s = [1,2,3,4]; type(s)
    <... 'list'>
    sage: s = (1,2,3,4); type(s)
    <... 'tuple'>
    sage: s = int(2006); type(s)
    <... 'int'>
    sage: s = float(2006); type(s)
    <... 'float'>

除此之外，Sage 还添加了许多其他类型。例如，向量空间：

::

    sage: V = VectorSpace(QQ, 1000000); V
    Vector space of dimension 1000000 over Rational Field
    sage: type(V)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category'>

只有某些函数可以在 ``V`` 上调用。
在其他数学软件系统中，这些函数可以使用“函数”符号 ``foo(V,...)`` 来调用。
在 Sage 中，某些函数附加到 ``V`` 类型（或类），并使用与 Java 或 C++ 类似的面向对象语法调用，
例如 ``V.foo(...)``。这种方式有助于保持全局命名空间的整洁，并允许名称相同但行为不同的函数存在，
而无需通过参数类型检查（或 case 语句）来决定调用哪个函数。此外，如果你重复使用函数名，该函数仍然可用
（例如，如果你调用某个函数 ``zeta``，那么要计算 Riemann-Zeta 函数在 0.5 处的值，
可以输入 ``s=.5; s.zeta()``）。

::

    sage: zeta = -1
    sage: s=.5; s.zeta()
    -1.46035450880959

在某些非常常见的情况下，为了方便起见，
同时避免使用面向对象符号可能导致数学表达式看起来令人困惑，
Sage 也支持常规的函数符号。这里有一些例子。

::

    sage: n = 2; n.sqrt()
    sqrt(2)
    sage: sqrt(2)
    sqrt(2)
    sage: V = VectorSpace(QQ,2)
    sage: V.basis()
        [(1, 0), (0, 1)]
    sage: basis(V)
        [(1, 0), (0, 1)]
    sage: M = MatrixSpace(GF(7), 2); M
    Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
    sage: A = M([1,2,3,4]); A
    [1 2]
    [3 4]
    sage: A.charpoly('x')
    x^2 + 2*x + 5
    sage: charpoly(A, 'x')
    x^2 + 2*x + 5

要列出 :math:`A` 的所有成员函数，请使用 tab 补全功能。
只需输入 ``A.``，然后在键盘上按 ``[tab]`` 键即可，
如 :ref:`section-tabcompletion` 中所述。

列表、元组和序列
============================

列表数据类型具有存储任意类型元素的功能。
和 C、C++ 等语言类似（但与大多数标准的计算机代数系统不同），列表元素的索引是从 :math:`0` 开始的：

::

    sage: v = [2, 3, 5, 'x', SymmetricGroup(3)]; v
    [2, 3, 5, 'x', Symmetric group of order 3! as a permutation group]
    sage: type(v)
    <... 'list'>
    sage: v[0]
    2
    sage: v[2]
    5

（在列表中进行索引时，不一定需要使用 Python 的整型作为索引！）
Sage 整数（或有理数，或任何具有 ``__index__`` 方法的对象）都可以正常使用。

::

    sage: v = [1,2,3]
    sage: v[2]
    3
    sage: n = 2      # Sage Integer
    sage: v[n]       # Perfectly OK!
    3
    sage: v[int(n)]  # Also OK.
    3

``range`` 函数创建一个包含 Python 整型（而不是 Sage 整数）元素的列表：

::

    sage: list(range(1, 15))
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

该函数在使用列表推导式构造列表时非常有用：

::

    sage: L = [factor(n) for n in range(1, 15)]
    sage: L
    [1, 2, 3, 2^2, 5, 2 * 3, 7, 2^3, 3^2, 2 * 5, 11, 2^2 * 3, 13, 2 * 7]
    sage: L[12]
    13
    sage: type(L[12])
    <class 'sage.structure.factorization_integer.IntegerFactorization'>
    sage: [factor(n) for n in range(1, 15) if is_odd(n)]
    [1, 3, 5, 7, 3^2, 11, 13]

有关如何使用列表推导式创建列表的更多内容，请参考 [PyT]_ 。

列表切片是一个非常好的功能。如果 ``L`` 是一个列表，那么 ``L[m:n]`` 返回
从第 :math:`m` 个元素开始到第 :math:`(n-1)` 个元素结束的子列表，如下所示：

::

    sage: L = [factor(n) for n in range(1, 20)]
    sage: L[4:9]
    [5, 2 * 3, 7, 2^3, 3^2]
    sage: L[:4]
    [1, 2, 3, 2^2]
    sage: L[14:4]
    []
    sage: L[14:]
    [3 * 5, 2^4, 17, 2 * 3^2, 19]

元组与列表类似，只不过它是不可变的，一旦创建便不能更改。

::

    sage: v = (1,2,3,4); v
    (1, 2, 3, 4)
    sage: type(v)
    <... 'tuple'>
    sage: v[1] = 5
    Traceback (most recent call last):
    ...
    TypeError: 'tuple' object does not support item assignment

序列是第三种面向列表的 Sage 类型。与列表和元组不同，序列不是 Python 内置类型。
默认情况下，序列是可变的，但可以使用 ``Sequence`` 类方法 ``set_immutable`` 将其设置为不可变，
如以下例子所示。序列的所有元素都属于同一个父对象，称为序列的领域 (universe)。

::

    sage: v = Sequence([1,2,3,4/5])
    sage: v
    [1, 2, 3, 4/5]
    sage: type(v)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: type(v[1])
    <class 'sage.rings.rational.Rational'>
    sage: v.universe()
    Rational Field
    sage: v.is_immutable()
    False
    sage: v.set_immutable()
    sage: v[0] = 3
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

序列派生自列表，可以在任何需要列表的地方使用：

::

    sage: v = Sequence([1,2,3,4/5])
    sage: isinstance(v, list)
    True
    sage: list(v)
    [1, 2, 3, 4/5]
    sage: type(list(v))
    <... 'list'>

另一个例子是，向量空间的基是不可变序列，因为不能改变它们至关重要。

::

    sage: V = QQ^3; B = V.basis(); B
    [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    sage: type(B)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: B[0] = B[1]
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.
    sage: B.universe()
    Vector space of dimension 3 over Rational Field

字典
============

字典（有时也被称为关联数组）是从“可哈希”对象（例如字符串、数字和元组等；详情请参见 Python 文档
http://docs.python.org/tut/node7.html 和
http://docs.python.org/lib/typesmapping.html ）
到任意对象的映射。

::

    sage: d = {1:5, 'sage':17, ZZ:GF(7)}
    sage: type(d)
    <... 'dict'>
    sage: list(d.keys())
    [1, 'sage', Integer Ring]
    sage: d['sage']
    17
    sage: d[ZZ]
    Finite Field of size 7
    sage: d[1]
    5

第三个键说明字典的索引可以很复杂，例如整数环。

你可以将上述字典转换为具有相同数据的列表：

.. link

::

    sage: list(d.items())
    [(1, 5), ('sage', 17), (Integer Ring, Finite Field of size 7)]

一种常见用法是遍历字典中的键值对：

::

    sage: d = {2:4, 3:9, 4:16}
    sage: [a*b for a, b in d.items()]
    [8, 27, 64]

正如最后的输出所示，字典是无序的。

集合
====

Python 有内建的集合类型。它提供的主要功能是快速查找元素是否在集合中，以及标准集合论运算。

::

    sage: X = set([1,19,'a']);   Y = set([1,1,1, 2/3])
    sage: X   # random sort order
    {1, 19, 'a'}
    sage: X == set(['a', 1, 1, 19])
    True
    sage: Y
    {2/3, 1}
    sage: 'a' in X
    True
    sage: 'a' in Y
    False
    sage: X.intersection(Y)
    {1}

Sage 也有自己的集合类型（在某些情况下使用 Python 内建集合类型实现），
但具有一些与 Sage 相关的额外功能。使用 ``Set(...)`` 来创建 Sage 集合。例如：

::

    sage: X = Set([1,19,'a']);   Y = Set([1,1,1, 2/3])
    sage: X   # random sort order
    {'a', 1, 19}
    sage: X == Set(['a', 1, 1, 19])
    True
    sage: Y
    {1, 2/3}
    sage: X.intersection(Y)
    {1}
    sage: print(latex(Y))
    \left\{1, \frac{2}{3}\right\}
    sage: Set(ZZ)
    Set of elements of Integer Ring

迭代器
=========

迭代器是 Python 最近添加的功能，在数学应用中特别有用。
这里有几个例子；详情请参见 [PyT]_ 。我们创建一个非负整数平方的迭代器，上限为 :math:`10000000`。

::

    sage: v = (n^2 for n in range(10000000))
    sage: next(v)
    0
    sage: next(v)
    1
    sage: next(v)
    4

我们创建一个 :math:`4p+1` 形式的素数迭代器，其中 :math:`p` 也是素数，并查看前几个值。

::

    sage: w = (4*p + 1 for p in Primes() if is_prime(4*p+1))
    sage: w         # in the next line, 0xb0853d6c is a random 0x number
    <generator object at 0xb0853d6c>
    sage: next(w)
    13
    sage: next(w)
    29
    sage: next(w)
    53

某些环，例如有限域和整数环有与之关联的迭代器：

::

    sage: [x for x in GF(7)]
    [0, 1, 2, 3, 4, 5, 6]
    sage: W = ((x,y) for x in ZZ for y in ZZ)
    sage: next(W)
    (0, 0)
    sage: next(W)
    (0, 1)
    sage: next(W)
    (0, -1)

循环、函数、控制语句和比较
=====================================================

我们已经看过了一些常见的 ``for`` 循环用法示例。在 Python 中，``for`` 循环具有缩进结构，例如：

.. CODE-BLOCK:: pycon

    >>> for i in range(5):
    ...     print(i)
    ...
    0
    1
    2
    3
    4

请注意 for 语句末尾的冒号（不像 GAP 或 Maple 中有 "do" 或 "od"），
以及循环体（即 ``print(i)``）前的缩进。这个缩进非常重要。
在 Sage 中，当你在 ":" 后按下 ``enter`` 时，会自动添加缩进，如下所示。

::

    sage: for i in range(5):
    ....:     print(i)  # now hit enter twice
    ....:
    0
    1
    2
    3
    4


符号 ``=`` 用于赋值。
符号 ``==`` 用于检查相等：

::

    sage: for i in range(15):
    ....:     if gcd(i,15) == 1:
    ....:         print(i)
    ....:
    1
    2
    4
    7
    8
    11
    13
    14

请牢记缩进如何决定 ``if``, ``for`` 和 ``while`` 语句的块结构：

::

    sage: def legendre(a,p):
    ....:     is_sqr_modp=-1
    ....:     for i in range(p):
    ....:         if a % p == i^2 % p:
    ....:             is_sqr_modp=1
    ....:     return is_sqr_modp

    sage: legendre(2,7)
    1
    sage: legendre(3,7)
    -1

当然，这不是勒让德符号 (Legendre symbol) 的高效实现！
它只是为了说明 Python/Sage 编程的各个方面。Sage 附带的函数 {kronecker}，
可以通过调用 PARI 的 C 库高效地计算勒让德符号。

最后，我们注意到数字之间的比较，如 ``==``, ``!=``, ``<=``, ``>=``, ``>``, ``<``，
会自动将两个数字转换为相同类型（如果可能的话）：

::

    sage: 2 < 3.1; 3.1 <= 1
    True
    False
    sage: 2/3 < 3/2;   3/2 < 3/1
    True
    True

使用 bool 来判断符号不等式：

::

    sage: x < x + 1
    x < x + 1
    sage: bool(x < x + 1)
    True

在比较不同类型的对象时，在大多数情况下，Sage 会尝试找到两者的共同复结构
（参见 :ref:`section-coercion` 了解更多细节）。
如果成功，比较将在强制转换的对象之间进行；如果不成功，则认为对象不相等。
要测试两个变量是否引用同一个对象，请使用 ``is``。
在下面这个示例中我们将看到，Python 整型 ``1`` 是唯一的，而 Sage 整型 ``1`` 则不是：

::

    sage: 1 is 2/2
    False
    sage: 1 is 1
    False
    sage: 1 == 2/2
    True

在以下两行代码中，第一个等式为 ``False``，因为没有从 :math:`\QQ \to \GF{5}` 的标准同态，
因此无法将 :math:`\GF{5}` 中的 :math:`1` 与 :math:`1 \in \QQ` 进行比较。
相反，由于存在从 :math:`\ZZ \to \GF{5}` 的标准映射，因此第二个比较为 ``True``。
需要注意的是，顺序不影响结果。

::

    sage: GF(5)(1) == QQ(1); QQ(1) == GF(5)(1)
    False
    False
    sage: GF(5)(1) == ZZ(1); ZZ(1) == GF(5)(1)
    True
    True
    sage: ZZ(1) == QQ(1)
    True

警告: Sage 中的比较比 Magma 更严格，Magma 会声明 :math:`1 \in \GF{5}` 等于 :math:`1 \in \QQ`。

::

    sage: magma('GF(5)!1 eq Rationals()!1')            # optional - magma
    true

性能分析
=========

    “过早优化乃万恶之源。” - Donald Knuth

.. sectionauthor:: Martin Albrecht <malb@informatik.uni-bremen.de>

有时检查代码中的瓶颈有助于了解哪些部分占用最多的计算时间；
这可以很好地了解哪些部分需要优化。
Python（以及 Sage）提供了几种性能分析工具和方法，
这个过程称为性能分析。

最简单的方式是使用交互式 shell 中的 ``prun`` 命令。
它会返回一个总结，描述哪些函数花了多少计算时间。
例如，要分析有限域上的矩阵乘法（版本 1.0 当前很慢！），可以这样做：

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])

.. skip

::

    sage: %prun B = A*A
           32893 function calls in 1.100 CPU seconds

    Ordered by: internal time

    ncalls tottime percall cumtime percall filename:lineno(function)
     12127  0.160   0.000   0.160  0.000 :0(isinstance)
      2000  0.150   0.000   0.280  0.000 matrix.py:2235(__getitem__)
      1000  0.120   0.000   0.370  0.000 finite_field_element.py:392(__mul__)
      1903  0.120   0.000   0.200  0.000 finite_field_element.py:47(__init__)
      1900  0.090   0.000   0.220  0.000 finite_field_element.py:376(__compat)
       900  0.080   0.000   0.260  0.000 finite_field_element.py:380(__add__)
         1  0.070   0.070   1.100  1.100 matrix.py:864(__mul__)
      2105  0.070   0.000   0.070  0.000 matrix.py:282(ncols)
      ...

这里 ``ncalls`` 是调用次数，``tottime`` 是给定函数花费的总时间（不包括调用子函数的时间），
``percall`` 是 ``tottime`` 除以 ``ncalls`` 的商。
``cumtime`` 是该函数及所有子函数花费的总时间（即，从调用到退出），
``percall`` 是 ``cumtime`` 除以原始调用次数的商，
``filename:lineno(function)`` 提供了每个函数的相关数据。
性能分析中的经验法则是：列表越靠前的函数，其代价越高，因而更需要进行优化。

与以往一样，``prun?`` 命令提供了使用性能分析器和理解输出详细信息的帮助。

性能分析数据还可以保存到一个对象中，以便进行更详细的检查：

.. skip

::

    sage: %prun -r A*A
    sage: stats = _
    sage: stats?

注意：输入 ``stats = prun -r A\*A`` 会显示语法错误消息，
因为 prun 是 IPython shell 命令而不是常规函数。

为了更好地以图形化方式呈现分析数据，你可以使用 hotshot 分析器，
``hotshot2cachetree`` 脚本，以及 ``kcachegrind`` 程序（仅限 Unix）。
以下是使用 hotshot 分析器的示例：

.. skip

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])
    sage: import hotshot
    sage: filename = "pythongrind.prof"
    sage: prof = hotshot.Profile(filename, lineevents=1)

.. skip

::

    sage: prof.run("A*A")
    <hotshot.Profile instance at 0x414c11ec>
    sage: prof.close()

这会在当前工作目录中生成一个 ``pythongrind.prof`` 文件。
现在可以将其转换为 cachegrind 格式进行可视化展示。

在系统终端中，输入

.. CODE-BLOCK:: shell-session

    $ hotshot2calltree -o cachegrind.out.42 pythongrind.prof

现在，输出文件 ``cachegrind.out.42`` 可以用 ``kcachegrind`` 查看。
请注意，需要遵守命名约定 ``cachegrind.out.XX``。
