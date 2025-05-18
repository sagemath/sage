.. _section-poly:

多项式
===========

在本节中，我们将介绍如何在 Sage 中创建和使用多项式。


.. _section-univariate:

一元多项式
----------------------

创建多项式环有三种方法。

::

    sage: R = PolynomialRing(QQ, 't')
    sage: R
    Univariate Polynomial Ring in t over Rational Field

这会创建一个多项式环，并告诉 Sage 在显示时使用字符串 't' 作为不定元。
然而，这并没有定义符号 ``t``，因此你不能用它来输入属于 ``R`` 的多项式（例如 :math:`t^2+1`）。

另一种方法是

.. link

::

    sage: S = QQ['t']
    sage: S == R
    True

这样做对于 ``t`` 也存在同样的问题。

第三种非常方便的方法是

::

    sage: R.<t> = PolynomialRing(QQ)

或

::

    sage: R.<t> = QQ['t']

甚至

::

    sage: R.<t> = QQ[]

这样做还有一个额外的好处，即它定义了变量 ``t`` 作为多项式环的不定元，
因此你可以轻松地构造 ``R`` 的元素，如下所示。
（请注意，第三种方法与 Magma 中的构造符号非常相似，并且可以像在 Magma 中一样用于广泛的对象。）

.. link

::

    sage: poly = (t+1) * (t+2); poly
    t^2 + 3*t + 2
    sage: poly in R
    True

无论你使用哪种方法定义多项式环，你都可以通过 :math:`0^{th}` 生成器恢复不定元：

::

    sage: R = PolynomialRing(QQ, 't')
    sage: t = R.0
    sage: t in R
    True

请注意，类似的构造方法适用于复数：复数可以被视为由符号 ``i`` 在实数上生成的，因此我们有以下内容：

::

    sage: CC
    Complex Field with 53 bits of precision
    sage: CC.0  # 0th generator of CC
    1.00000000000000*I

对于多项式环，你可以在创建环时同时获得环及其生成器，或者仅获得生成器，如下所示：

::

    sage: R, t = QQ['t'].objgen()
    sage: t    = QQ['t'].gen()
    sage: R, t = objgen(QQ['t'])
    sage: t    = gen(QQ['t'])

最后我们在 :math:`\QQ[t]` 中进行一些算术运算。

::

    sage: R, t = QQ['t'].objgen()
    sage: f = 2*t^7 + 3*t^2 - 15/19
    sage: f^2
    4*t^14 + 12*t^9 - 60/19*t^7 + 9*t^4 - 90/19*t^2 + 225/361
    sage: cyclo = R.cyclotomic_polynomial(7); cyclo
    t^6 + t^5 + t^4 + t^3 + t^2 + t + 1
    sage: g = 7 * cyclo * t^5 * (t^5 + 10*t + 2)
    sage: g
    7*t^16 + 7*t^15 + 7*t^14 + 7*t^13 + 77*t^12 + 91*t^11 + 91*t^10 + 84*t^9
           + 84*t^8 + 84*t^7 + 84*t^6 + 14*t^5
    sage: F = factor(g); F
    (7) * t^5 * (t^5 + 10*t + 2) * (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1)
    sage: F.unit()
    7
    sage: list(F)
    [(t, 5), (t^5 + 10*t + 2, 1), (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1, 1)]

注意，因式分解正确考虑并记录了单位部分。

如果你在某个研究项目中大量使用某个函数，例如 ``R.cyclotomic_polynomial``，
除了引用 Sage 之外，你还应该尝试找出 Sage 的哪个组件在实际计算分圆多项式并引用它。
在这种情况下，如果你输入 ``R.cyclotomic_polynomial??`` 查看源代码，
你很快会看到一行 ``f = pari.polcyclo(n)``，这意味着 PARI 被用于计算分圆多项式。
你的作品中也需要引用 PARI。

除以两个多项式会构造分数域的元素（Sage 会自动创建）。

::

    sage: x = QQ['x'].0
    sage: f = x^3 + 1; g = x^2 - 17
    sage: h = f/g;  h
    (x^3 + 1)/(x^2 - 17)
    sage: h.parent()
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

使用 Laurent 级数，可以在 ``QQ[x]`` 的分数域中计算级数展开：

::

    sage: R.<x> = LaurentSeriesRing(QQ); R
    Laurent Series Ring in x over Rational Field
    sage: 1/(1-x) + O(x^10)
    1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)

如果我们给变量不同的命名，我们会得到不同的一元多项式环。

::

    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<y> = PolynomialRing(QQ)
    sage: x == y
    False
    sage: R == S
    False
    sage: R(y)
    x
    sage: R(y^2 - 17)
    x^2 - 17

环由变量决定。请注意，使用名为 ``x`` 的变量创建另一个环不会返回不同的环。

::

    sage: R = PolynomialRing(QQ, "x")
    sage: T = PolynomialRing(QQ, "x")
    sage: R == T
    True
    sage: R is T
    True
    sage: R.0 == T.0
    True

Sage 还支持任意基环上的幂级数和 Laurent 级数环。
在下面的示例中，我们创建了 `\GF{7}[[T]]` 的一个元素，
并通过相除创建 :math:`\GF{7}((T))` 的一个元素。

::

    sage: R.<T> = PowerSeriesRing(GF(7)); R
    Power Series Ring in T over Finite Field of size 7
    sage: f = T  + 3*T^2 + T^3 + O(T^4)
    sage: f^3
    T^3 + 2*T^4 + 2*T^5 + O(T^6)
    sage: 1/f
    T^-1 + 4 + T + O(T^2)
    sage: parent(1/f)
    Laurent Series Ring in T over Finite Field of size 7

你也可以使用双括号简写来创建幂级数环：

::

    sage: GF(7)[['T']]
    Power Series Ring in T over Finite Field of size 7

多元多项式
------------------------

要处理多个变量的多项式，我们首先声明多项式环和变量。

::

    sage: R = PolynomialRing(GF(5),3,"z") # here, 3 = number of variables
    sage: R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

与定义一元多项式环一样，有多种方法：

::

    sage: GF(5)['z0, z1, z2']
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5
    sage: R.<z0,z1,z2> = GF(5)[]; R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

此外，如果你想让变量名为单个字母，你可以使用以下简写：

::

    sage: PolynomialRing(GF(5), 3, 'xyz')
    Multivariate Polynomial Ring in x, y, z over Finite Field of size 5

接下来让我们进行一些算术运算。

::

    sage: z = GF(5)['z0, z1, z2'].gens()
    sage: z
    (z0, z1, z2)
    sage: (z[0]+z[1]+z[2])^2
    z0^2 + 2*z0*z1 + z1^2 + 2*z0*z2 + 2*z1*z2 + z2^2

你还可以使用更多数学符号来构造多项式环。

::

    sage: R = GF(5)['x,y,z']
    sage: x,y,z = R.gens()
    sage: QQ['x']
    Univariate Polynomial Ring in x over Rational Field
    sage: QQ['x,y'].gens()
    (x, y)
    sage: QQ['x'].objgens()
    (Univariate Polynomial Ring in x over Rational Field, (x,))

多元多项式在 Sage 中使用 Python 字典和多项式的“分配表示”实现。
Sage 使用了一些 Singular [Si]_ ，例如，用于计算理想的最大公约数和 Gröbner 基。

::

    sage: R, (x, y) = PolynomialRing(RationalField(), 2, 'xy').objgens()
    sage: f = (x^3 + 2*y^2*x)^2
    sage: g = x^2*y^2
    sage: f.gcd(g)
    x^2

接下来我们通过简单地将 ``(f,g)`` 乘以 ``R``
来创建由 :math:`f` 和 :math:`g` 生成的理想 :math:`(f,g)`，（也可以写做 ``ideal([f,g])`` 或 ``ideal(f,g)``）。

.. link

::

    sage: I = (f, g)*R; I
    Ideal (x^6 + 4*x^4*y^2 + 4*x^2*y^4, x^2*y^2) of Multivariate Polynomial
    Ring in x, y over Rational Field
    sage: B = I.groebner_basis(); B
    [x^6, x^2*y^2]
    sage: x^2 in I
    False

顺便说一句，上面的 Gröbner 基不是一个列表，而是一个不可变序列。
这意味着它有全集，父结构，并且不可更改（这是好的，因为更改基会破坏使用 Gröbner 基的其他例程）。

.. link

::

    sage: B.universe()
    Multivariate Polynomial Ring in x, y over Rational Field
    sage: B[1] = x
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Sage 中有一些（没有我们想要的那么多）交换代数可用，通过 Singular 实现。
例如，我们可以计算 :math:`I` 的初等分解和相关素数：

.. link

::

    sage: I.primary_decomposition()
    [Ideal (x^2) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y^2, x^6) of Multivariate Polynomial Ring in x, y over Rational Field]
    sage: I.associated_primes()
    [Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field]
