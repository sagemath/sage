
赋值、等式和算术
====================================

Sage 基本上使用 Python 编程语言，因此大多数 Python 入门书籍都能帮助你学习 Sage。

Sage 使用 ``=`` 进行赋值。使用 ``==``, ``<=``, ``>=``, ``<`` 和 ``>`` 进行比较：

::

    sage: a = 5
    sage: a
    5
    sage: 2 == 2
    True
    sage: 2 == 3
    False
    sage: 2 < 3
    True
    sage: a == 5
    True

Sage 提供所有基本的数学运算：

::

    sage: 2**3    #  ** means exponent
    8
    sage: 2^3     #  ^ is a synonym for ** (unlike in Python)
    8
    sage: 10 % 3  #  for integer arguments, % means mod, i.e., remainder
    1
    sage: 10/4
    5/2
    sage: 10//4   #  for integer arguments, // returns the integer quotient
    2
    sage: 4 * (10 // 4) + 10 % 4 == 10
    True
    sage: 3^2*4 + 2%5
    38

像 ``3^2*4 + 2%5`` 这样的表达式的计算取决于运算的顺序；
:ref:`section-precedence` 中的“运算符优先级表”给出了明确的规定。

Sage 还提供了许多常见数学函数；以下是一些例子：

::

    sage: sqrt(3.4)
    1.84390889145858
    sage: sin(5.135)
    -0.912021158525540
    sage: sin(pi/3)
    1/2*sqrt(3)

如最后一个例子所示，一些数学表达式返回“精确”值，而不是近似值。
要获得数值近似，可以使用函数 ``N`` 或方法 ``n`` （二者都有一个更长的名称 ``numerical_approx``，函数 ``N`` 与 ``n`` 相同）。
这些函数接受可选参数 ``prec``，即请求的精度位数，以及 ``digits``，即请求的十进制精度位数；默认精度为 53 位。

::

    sage: exp(2)
    e^2
    sage: n(exp(2))
    7.38905609893065
    sage: sqrt(pi).numerical_approx()
    1.77245385090552
    sage: sin(10).n(digits=5)
    -0.54402
    sage: N(sin(10),digits=10)
    -0.5440211109
    sage: numerical_approx(pi, prec=200)
    3.1415926535897932384626433832795028841971693993751058209749

Python 是动态类型语言，所以每个变量引用的值都有一个类型与之关联，
但在给定作用域内，一个给定变量可以保存任意 Python 类型的值：

::

    sage: a = 5   # a is an integer
    sage: type(a)
    <class 'sage.rings.integer.Integer'>
    sage: a = 5/3  # now a is a rational number
    sage: type(a)
    <class 'sage.rings.rational.Rational'>
    sage: a = 'hello'  # now a is a string
    sage: type(a)
    <... 'str'>

C 语言作为静态类型语言就非常不同；一个声明为 int 类型的变量在其作用域内只能保存 int。
