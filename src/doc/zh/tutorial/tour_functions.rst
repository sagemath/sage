.. _section-functions-issues:

常见函数问题
=================================

定义函数的某些方面（例如，用于微分或绘图）可能会令人困惑。我们尝试在本节中解答一些相关问题。

以下是几种可以被称为“函数”的定义方法：

1. 定义一个 Python 函数，如 :ref:`section-functions` 中所述。
这些函数可以被绘制，但不能被微分或积分。

::

       sage: def f(z): return z^2
       sage: type(f)
       <... 'function'>
       sage: f(3)
       9
       sage: plot(f, 0, 2)
       Graphics object consisting of 1 graphics primitive

请注意最后一行的语法。使用 ``plot(f(z), 0, 2)`` 会报 :class:`NameError`。
因为 ``z`` 是 ``f`` 定义中的一个虚拟变量，在该定义之外未定义。
为了能够在 plot 命令中使用 ``f(z)``，需要将 ``z`` （或其他所需内容）定义为变量。
我们可以使用下面的语法，或者采用我们给出的第二种方法。

.. link

::

       sage: var('z')   # define z to be a variable
       z
       sage: f(z)
       z^2
       sage: plot(f(z), 0, 2)
       Graphics object consisting of 1 graphics primitive

此时，``f(z)`` 是一个符号表达式，即我们接下来要介绍的方法。

2. 定义一个“可调用的符号表达式”。这些表达式可以被绘制、微分和积分。

::

       sage: g(x) = x^2
       sage: g        # g sends x to x^2
       x |--> x^2
       sage: g(3)
       9
       sage: Dg = g.derivative(); Dg
       x |--> 2*x
       sage: Dg(3)
       6
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(g, 0, 2)
       Graphics object consisting of 1 graphics primitive

注意，虽然 ``g`` 是一个可调用的符号表达式，但 ``g(x)`` 是一个相关但不同类型的对象，
尽管存在一些问题，但它也可以被绘制、微分等：请参见下文中的第 5 点。

.. link

::

       sage: g(x)
       x^2
       sage: type(g(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: g(x).derivative()
       2*x
       sage: plot(g(x), 0, 2)
       Graphics object consisting of 1 graphics primitive

3. 使用预定义的 Sage '微积分函数'。 这些函数可以被绘制，并且稍加辅助可以进行微分和积分。

::

       sage: type(sin)
       <class 'sage.functions.trig.Function_sin'>
       sage: plot(sin, 0, 2)
       Graphics object consisting of 1 graphics primitive
       sage: type(sin(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(sin(x), 0, 2)
       Graphics object consisting of 1 graphics primitive

单独使用 ``sin`` 不能被微分，至少不能得到 ``cos``.

::

       sage: f = sin
       sage: f.derivative()
       Traceback (most recent call last):
       ...
       AttributeError: ...

用 ``f = sin(x)`` 替换 ``sin`` 就可以了，
但更好的方法可能是使用 ``f(x) = sin(x)`` 来定义一个可调用的符号表达式。

::

       sage: S(x) = sin(x)
       sage: S.derivative()
       x |--> cos(x)

以下是一些常见问题及其解释：

\4. 非预期执行

::

       sage: def h(x):
       ....:     if x<2:
       ....:         return 0
       ....:     else:
       ....:         return x-2


问题： ``plot(h(x), 0, 4)`` 绘制的是直线 `y=x-2`，而不是由 ``h`` 定义的分段函数。
原因是，在命令 ``plot(h(x), 0, 4)`` 中，首先执行 ``h(x)``，这意味着将符号变量 ``x`` 插入函数 ``h`` 中。
因此，不等式 ``x < 2`` 首先执行得到 ``False``，因此 ``h(x)`` 会执行 ``x - 2``。
可以通过以下方法看到这个过程

.. link

::

        sage: bool(x < 2)
        False
        sage: h(x)
        x - 2

注意，这里有两个不同的 ``x``：用于定义函数 ``h`` 的 Python 变量（在其定义中是局部的）和 Sage 启动时可用的符号变量 ``x``。

解决方案：不要使用 ``plot(h(x), 0, 4)``；而是使用

.. link

::

       sage: plot(h, 0, 4)
       Graphics object consisting of 1 graphics primitive

\5. 意外产生常数而非函数。

::

       sage: f = x
       sage: g = f.derivative()
       sage: g
       1

问题：以 ``g(3)`` 为例，会返回一个错误，
提示 "ValueError: the number of arguments must be less than or equal to 0."。

.. link

::

       sage: type(f)
       <class 'sage.symbolic.expression.Expression'>
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>

``g`` 不是函数，而是一个常数，所以它没有关联的变量，不能将任何内容插入其中。

解决方案：有几种选择。

- 将 ``f`` 定义为符号表达式。

::

         sage: f(x) = x        # instead of 'f = x'
         sage: g = f.derivative()
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>

- 或者保留 ``f`` 的原始定义，将 ``g`` 定义为符号表达式。

::

         sage: f = x
         sage: g(x) = f.derivative()  # instead of 'g = f.derivative()'
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>

- 抑或保留 ``f`` 和 ``g`` 的原始定义，指定要替换的变量。

::

         sage: f = x
         sage: g = f.derivative()
         sage: g
         1
         sage: g(x=3)    # instead of 'g(3)'
         1

最后，还有另一种方法可以区分 ``f = x`` 和 ``f(x) = x`` 的导数

::

       sage: f(x) = x
       sage: g = f.derivative()
       sage: g.variables()  # the variables present in g
       ()
       sage: g.arguments()  # the arguments which can be plugged into g
       (x,)
       sage: f = x
       sage: h = f.derivative()
       sage: h.variables()
       ()
       sage: h.arguments()
       ()

正如上面例子试图说明的那样，``h`` 不接受任何参数，这就是为什么 ``h(3)`` 会返回错误。
