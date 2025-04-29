.. _section-rings:

基本环
===========

在定义矩阵、向量或多项式时，指定它们所定义的“环”非常有用，有时甚至是必须的。
*环* 是一种数学结构，具有良好的加法和乘法概念；如果你以前从未听说过它们，你可能只需要了解以下四种常用的环：

* 整数 `\{..., -1, 0, 1, 2, ...\}`，在 Sage 中称为 ``ZZ``。
* 有理数 -- 即分数或整数的比率 -- 在 Sage 中称为 ``QQ``。
* 实数，在 Sage 中称为 ``RR``。
* 复数，在 Sage 中称为 ``CC``。

了解这些区别是必要的，因为同一个多项式可能会根据它所定义的环而有所不同。
例如，多项式 `x^2-2` 有两个根，`\pm \sqrt{2}`。
这些根不是有理数，所以如果你处理的是具有有理系数的多项式，那么这个多项式无法因式分解。
但使用实系数，它便可以因式分解。
因此，你可能需要指定环以确保获得预期的信息。
以下两个命令分别定义了具有有理系数和实系数的多项式集。
集合被命名为 "ratpoly" 和 "realpoly"，但这里并不重要；
然而，请注意字符串 ".<t>" 和 ".<z>" 分别命名了两种情况下使用的 *变量*。::

    sage: ratpoly.<t> = PolynomialRing(QQ)
    sage: realpoly.<z> = PolynomialRing(RR)

现在我们来演示 `x^2-2` 的因式分解:

.. link

::

    sage: factor(t^2-2)
    t^2 - 2
    sage: factor(z^2-2)
    (z - 1.41421356237310) * (z + 1.41421356237310)

类似的情况也适用于矩阵：矩阵的行简化形式可能取决于它所定义的环，以及它的特征值和特征向量。
有关构造多项式的更多信息，请参见 :ref:`section-poly`，
有关矩阵的更多信息，请参见 :ref:`section-linalg`。

符号 ``I`` 表示 :math:`-1` 的平方根；``i`` 是 ``I`` 的同义词。显然，它不是一个有理数::

    sage: i  # square root of -1
    I
    sage: i in QQ
    False

注意：如果变量 ``i`` 已被赋予其他值，例如，如果它被用作循环变量，则上述代码可能无法按预期工作。如果是这种情况，请输入::

    sage: reset('i')

以获得 ``i`` 的原始复数值。

定义复数时有一个需要注意的地方：如上所述，符号 ``i`` 表示 `-1` 的平方根，
但是它是 `-1` 的*形式*平方根，是一个代数数。
调用 ``CC(i)`` 或 ``CC.0`` 或 ``CC.gen(0)`` 返回 `-1` 的*复数*平方根。
通过所谓的强制转换，可以进行涉及不同类型数字的算术运算，请参见 :ref:`section-coercion`。

::

    sage: i = CC(i)       # floating point complex number
    sage: i == CC.0
    True
    sage: a, b = 4/3, 2/3
    sage: z = a + b*i
    sage: z
    1.33333333333333 + 0.666666666666667*I
    sage: z.imag()        # imaginary part
    0.666666666666667
    sage: z.real() == a   # automatic coercion before comparison
    True
    sage: a + b
    2
    sage: 2*b == a
    True
    sage: parent(2/3)
    Rational Field
    sage: parent(4/2)
    Rational Field
    sage: 2/3 + 0.1       # automatic coercion before addition
    0.766666666666667
    sage: 0.1 + 2/3       # coercion rules are symmetric in Sage
    0.766666666666667

以下是 Sage 中一些基本环的更多示例。
如上所述，有理数环可以使用 ``QQ`` 或 ``RationalField()`` 来引用
（*域* 是满足乘法交换律的环，且每个非零元素在该环中都有一个倒数，因此有理数构成一个域，但整数不构成）::

    sage: RationalField()
    Rational Field
    sage: QQ
    Rational Field
    sage: 1/2 in QQ
    True

十进制数 ``1.2`` 被认为是 `QQ`` 中的数：
也可以“强制转换”成有理数的十进制数被认为是有理数（参见 :ref:`section-coercion`）。
数字 `\pi` 和 `\sqrt{2}` 不是有理数::

    sage: 1.2 in QQ
    True
    sage: pi in QQ
    False
    sage: pi in RR
    True
    sage: sqrt(2) in QQ
    False
    sage: sqrt(2) in CC
    True

为了在高等数学中使用，Sage 还具备其他环，例如有限域，`p`-adic 整数，代数数环，多项式环和矩阵环。
以下是其中一些的构造::

    sage: GF(3)
    Finite Field of size 3
    sage: GF(27, 'a')  # need to name the generator if not a prime field
    Finite Field in a of size 3^3
    sage: Zp(5)
    5-adic Ring with capped relative precision 20
    sage: sqrt(3) in QQbar # algebraic closure of QQ
    True
