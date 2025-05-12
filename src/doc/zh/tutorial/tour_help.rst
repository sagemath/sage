.. _chapter-help:

获取帮助
============

Sage 有大量的内置文档，可以通过输入函数或常量的名称，然后加上问号来访问：

.. skip

::

    sage: tan?
    Type:        <class 'sage.calculus.calculus.Function_tan'>
    Definition:  tan( [noargspec] )
    Docstring:

        The tangent function

        EXAMPLES:
            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
    sage: log2?
    Type:        <class 'sage.functions.constants.Log2'>
    Definition:  log2( [noargspec] )
    Docstring:

        The natural logarithm of the real number 2.

        EXAMPLES:
            sage: log2
            log2
            sage: float(log2)
            0.69314718055994529
            sage: RR(log2)
            0.693147180559945
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(log2)
            0.69314718055994530941723212145817656807550013436025525412068
            sage: l = (1-log2)/(1+log2); l
            (1 - log(2))/(log(2) + 1)
            sage: R(l)
            0.18123221829928249948761381864650311423330609774776013488056
            sage: maxima(log2)
            log(2)
            sage: maxima(log2).float()
            .6931471805599453
            sage: gp(log2)
            0.6931471805599453094172321215             # 32-bit
            0.69314718055994530941723212145817656807   # 64-bit
    sage: sudoku?
    File:        sage/local/lib/python2.5/site-packages/sage/games/sudoku.py
    Type:        <... 'function'>
    Definition:  sudoku(A)
    Docstring:

        Solve the 9x9 Sudoku puzzle defined by the matrix A.

        EXAMPLE:
            sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0,
        0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0, 0,0,0, 2,0,8, 0,0,0,
        0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2,
        0,0,0, 4,9,0, 0,5,0, 0,0,3])
            sage: A
            [5 0 0 0 8 0 0 4 9]
            [0 0 0 5 0 0 0 3 0]
            [0 6 7 3 0 0 0 0 1]
            [1 5 0 0 0 0 0 0 0]
            [0 0 0 2 0 8 0 0 0]
            [0 0 0 0 0 0 0 1 8]
            [7 0 0 0 0 4 1 5 0]
            [0 3 0 0 0 2 0 0 0]
            [4 9 0 0 5 0 0 0 3]
            sage: sudoku(A)
            [5 1 3 6 8 7 2 4 9]
            [8 4 9 5 2 1 6 3 7]
            [2 6 7 3 4 9 5 8 1]
            [1 5 8 4 6 3 9 7 2]
            [9 7 4 2 1 8 3 6 5]
            [3 2 6 7 9 5 4 1 8]
            [7 8 2 9 3 4 1 5 6]
            [6 3 5 1 7 2 8 9 4]
            [4 9 1 8 5 6 7 2 3]

Sage 还提供了“Tab 补全”功能：输入函数的前几个字母，然后按下 :kbd:`Tab` 键。
例如，如果你输入 ``ta`` 然后按下 :kbd:`Tab`，Sage 会显示 ``tachyon, tan, tanh,
taylor``。这是查找 Sage 中函数和其他结构名称的好方法。


.. _section-functions:

函数、缩进和计数
====================================

在 Sage 中定义一个新函数，请使用 ``def`` 命令，并在变量名列表后加上冒号。例如：

::

    sage: def is_even(n):
    ....:     return n%2 == 0
    sage: is_even(2)
    True
    sage: is_even(3)
    False

注意：根据你查看的教程版本，你可能会在本例的第二行看到三个点 ``....:``。
请勿输入它们，它们只是为了强调代码的缩进。
在这种情况下，请在块末尾按 [Return/Enter] 以插入空行并结束函数定义。

你不需要指定输入参数的类型。
你可以指定多个输入，每个输入都可以有一个可选的默认值。
例如，如果未指定 ``divisor``，则下面的函数默认值为 ``divisor=2``。

::

    sage: def is_divisible_by(number, divisor=2):
    ....:     return number%divisor == 0
    sage: is_divisible_by(6,2)
    True
    sage: is_divisible_by(6)
    True
    sage: is_divisible_by(6, 5)
    False

调用函数时，你还可以显式地指定一个或多个输入；如果你显式地指定输入，可以以任意顺序给出它们：

.. link

::

    sage: is_divisible_by(6, divisor=5)
    False
    sage: is_divisible_by(divisor=2, number=6)
    True

在 Python 中，代码块不是用大括号或其他语言中的开始和结束标记来表示的。
相反，代码块由缩进来表示，缩进必须完全匹配。
例如，以下是一个语法错误，因为 ``return`` 语句的缩进与上面的其他行不一致：

.. skip

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:    return v
    Syntax Error:
           return v

如果你修复了缩进，函数就可以正常工作：

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:     return v
    sage: even(10)
    [4, 6, 8]

行末不需要分号；在大多数情况下，行以换行符结束。但是，你可以在一行上放置多个语句，用分号间隔：

::

    sage: a = 5; b = a + 3; c = b^2; c
    64

如果你希望一行代码跨越多行，可以使用反斜杠：

::

    sage: 2 + \
    ....:    3
    5

在 Sage 中，你可以通过遍历整数区间来计数。
例如，下面代码的第一行与 C++ 或 Java 中的 ``for(i=0; i<3; i++)`` 完全一样：

::

    sage: for i in range(3):
    ....:     print(i)
    0
    1
    2

下面代码的第一行与 ``for(i=2;i<5;i++)`` 等价。

::

    sage: for i in range(2,5):
    ....:     print(i)
    2
    3
    4

第三个参数控制步长，所以下面代码与 ``for(i=1;i<6;i+=2)`` 等价。

::

    sage: for i in range(1,6,2):
    ....:     print(i)
    1
    3
    5

通常你会希望创建一个漂亮的表格来显示你使用 Sage 计算的数字。
一个简单的方法是使用格式化字符串。
下面，我们创建三个宽度正好为 6 的列，并制作一个平方和立方的表格。

::

    sage: for i in range(5):
    ....:     print('%6s %6s %6s' % (i, i^2, i^3))
         0      0      0
         1      1      1
         2      4      8
         3      9     27
         4     16     64

Sage 中最基本的数据结构是列表，顾名思义，就是一个任意对象的列表。
例如，以下命令使用 ``range`` 创建一个列表::

    sage: list(range(2,10))
    [2, 3, 4, 5, 6, 7, 8, 9]

下面是一个更复杂的列表：

::

    sage: v = [1, "hello", 2/3, sin(x^3)]
    sage: v
    [1, 'hello', 2/3, sin(x^3)]

如如许多编程语言一样，列表的索引是从 0 开始。

.. link

::

    sage: v[0]
    1
    sage: v[3]
    sin(x^3)

使用 ``len(v)`` 获取 ``v`` 的长度，
使用 ``v.append(obj)`` 将新对象追加到 ``v`` 的末尾，
使用 ``del v[i]`` 删除 ``v`` 的第 :math:`i` 项：

.. link

::

    sage: len(v)
    4
    sage: v.append(1.5)
    sage: v
    [1, 'hello', 2/3, sin(x^3), 1.50000000000000]
    sage: del v[1]
    sage: v
    [1, 2/3, sin(x^3), 1.50000000000000]

另一个重要的数据结构是字典（或关联数组）。
字典的工作方式类似于列表，但它可以用几乎任何对象来索引（索引必须是不可变的）：


::

    sage: d = {'hi':-2,  3/8:pi,   e:pi}
    sage: d['hi']
    -2
    sage: d[e]
    pi

你还可以使用类定义新的数据类型。
使用类封装数学对象是一种强大的技术，可以帮助简化和组织你的 Sage 程序。
下面，我们定义一个表示不超过 *n* 的正偶数列表的类；它从内置类型 ``list`` 派生而来。

::

    sage: class Evens(list):
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:         list.__init__(self, range(2, n+1, 2))
    ....:     def __repr__(self):
    ....:         return "Even positive numbers up to n."

``__init__`` 方法在创建对象时调用以初始化对象；
``__repr__`` 方法打印对象。
我们在 ``__init__`` 方法的第二行调用列表构造函数。
下面我们创建 ``Evens`` 类的对象：

.. link

::

    sage: e = Evens(10)
    sage: e
    Even positive numbers up to n.

注意，``e`` 使用我们定义的 ``__repr__`` 方法打印。
要查看底层数字列表，请使用 ``list`` 函数：

.. link

::

    sage: list(e)
    [2, 4, 6, 8, 10]

我们还可以访问属性 ``n`` 或像列表一样操作 ``e``。

.. link

::

    sage: e.n
    10
    sage: e[2]
    6
