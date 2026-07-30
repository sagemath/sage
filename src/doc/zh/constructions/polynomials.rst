******
多项式
******

.. index::
   pair: polynomial; powers

.. _section-polynomialpower:

多项式幂
========

如何在 Sage 中计算模多项式幂？

要计算 `GF(97)[x]` 中的 `x^{2006} \pmod {x^3 + 7}`，
我们需要创建商环 `GF(97)[x]/(x^3+7)`，
并在其中计算 `x^{2006}`。
作为 Sage 符号的一个重要细节，
我们必须区分 `GF(97)[x]` 中的 `x` 和商环 `GF(97)[x]/(x^3+7)` 中的的对应元素
（我们用 `a` 来表示）。

::

    sage: R = PolynomialRing(GF(97),'x')
    sage: x = R.gen()
    sage: S = R.quotient(x^3 + 7, 'a')
    sage: a = S.gen()
    sage: S
    Univariate Quotient Polynomial Ring in a over
    Finite Field of size 97 with modulus x^3 + 7
    sage: a^2006
    4*a^2

另一种计算方法是：

::

    sage: R = PolynomialRing(GF(97),'x')
    sage: x = R.gen()
    sage: S = R.quotient(x^3 + 7, 'a')
    sage: a = S.gen()
    sage: a^20062006
    80*a
    sage: libgap.eval("R:= PolynomialRing( GF(97))")
    GF(97)[x_1]
    sage: libgap.eval("i:= IndeterminatesOfPolynomialRing(R)")
    [ x_1 ]
    sage: libgap.eval("x:= i[1]"); libgap.eval("f:= x;")
    x_1
    x_1
    sage: libgap.eval("PowerMod( R, x, 20062006, x^3+7 );")
    Z(97)^41*x_1
    sage: libgap.eval("PowerMod( R, x, 20062006, x^3+7 );")
    Z(97)^41*x_1
    sage: libgap.eval("PowerMod( R, x, 2006200620062006, x^3+7 );")
    Z(97)^4*x_1^2
    sage: a^2006200620062006
    43*a^2
    sage: libgap.eval("PowerMod( R, x, 2006200620062006, x^3+7 );")
    Z(97)^4*x_1^2
    sage: libgap.eval("Int(Z(97)^4)")
    43

.. index::
   pair: polynomial; factorization

.. _section-factor:

因式分解
========

你可以使用 Sage 对多项式进行因式分解。

在 Sage 中对一元多项式进行因式分解，
只需将方法 ``factor`` 应用到 PolynomialRingElement 对象 f 上。
实际上，该方法会调用 Pari，因此计算相当快。

::

    sage: x = PolynomialRing(RationalField(), 'x').gen()
    sage: f = (x^3 - 1)^2-(x^2-1)^2
    sage: f.factor()
    (x - 1)^2 * x^2 * (x^2 + 2*x + 2)

使用 Singular 接口，Sage 还可以对多元多项式进行因式分解。

::

    sage: x, y = PolynomialRing(RationalField(), 2, ['x','y']).gens()
    sage: f =  (9*y^6 - 9*x^2*y^5 - 18*x^3*y^4 - 9*x^5*y^4 + 9*x^6*y^2 + 9*x^7*y^3
    ....:     + 18*x^8*y^2 - 9*x^11)
    sage: f.factor()
    (9) * (-x^5 + y^2) * (x^6 - 2*x^3*y^2 - x^2*y^3 + y^4)

.. index::
   pair: polynomial; gcd

多项式最大公约数
================

下面这个例子展示了一元多项式最大公约数：

::

    sage: x = PolynomialRing(RationalField(), 'x').gen()
    sage: f = 3*x^3 + x
    sage: g = 9*x*(x+1)
    sage: f.gcd(g)
    x

下面这个例子展示了多元多项式最大公约数：

::

    sage: R.<x,y,z> = PolynomialRing(RationalField(), order='lex')
    sage: f = 3*x^2*(x+y)
    sage: g = 9*x*(y^2 - x^2)
    sage: f.gcd(g)
    x^2 + x*y

下面是另一种方法：

::

    sage: R2 = singular.ring(0, '(x,y,z)', 'lp')
    sage: a = singular.new('3x2*(x+y)')
    sage: b = singular.new('9x*(y2-x2)')
    sage: g = a.gcd(b)
    sage: g
    x^2+x*y

下面这个例子展示了通过 GAP 接口计算一元多项式最大公约数。

::

    sage: R = libgap.PolynomialRing(GF(2)); R
    GF(2)[x_1]
    sage: i = R.IndeterminatesOfPolynomialRing(); i
    [ x_1 ]
    sage: x_1 = i[0]
    sage: f = (x_1^3 - x_1 + 1)*(x_1 + x_1^2); f
    x_1^5+x_1^4+x_1^3+x_1
    sage: g = (x_1^3 - x_1 + 1)*(x_1 + 1); g
    x_1^4+x_1^3+x_1^2+Z(2)^0
    sage: f.Gcd(g)
    x_1^4+x_1^3+x_1^2+Z(2)^0

当然，我们也可以在生成器上执行相同的计算，
它使用 NTL 库（该库能够非常快速地处理有限域上的大规模多项式最大公约数计算）。

::

    sage: x = PolynomialRing(GF(2), 'x').gen()
    sage: f = (x^3 - x + 1)*(x + x^2); f
    x^5 + x^4 + x^3 + x
    sage: g = (x^3 - x + 1)*(x + 1)
    sage: f.gcd(g)
    x^4 + x^3 + x^2 + 1

.. index::
   pair: polynomial; roots

.. _section-roots:

多项式的根
==========

Sage 可以计算一元多项式的根。

::

    sage: x = PolynomialRing(RationalField(), 'x').gen()
    sage: f = x^3 - 1
    sage: f.roots()
    [(1, 1)]
    sage: f = (x^3 - 1)^2
    sage: f.roots()
    [(1, 2)]
    sage: x = PolynomialRing(CyclotomicField(3), 'x').gen()
    sage: f = x^3 - 1
    sage: f.roots()
    [(1, 1), (zeta3, 1), (-zeta3 - 1, 1)]

第一个元素是根，第二个元素是它的重数。

在某些情况下，GAP 确实可以求解一元多项式的根，
但 GAP 通常不会这样做（根必须生成有限域或循环域的子域）。
然而，有一个名为 ``RadiRoot`` 的 GAP 包，必须将其安装到 GAP 中，
因为它确实有助于为有理系数的多项式执行此操作（``radiroot`` 本身需要安装其他包；请参阅其网页了解更多详情）。
``Factors`` 命令实际上有一个选项，允许你增加基域，以便因式分解实际返回根。
更多详情，请参阅 GAP 参考手册第 64.10 节“多项式因式分解”中给出的示例。

.. index::
   pair: polynomial; evaluation

.. _section-evaluate:

多元函数求值
============

你可以像往常一样在 Sage 中通过代入点来计算多项式的值：

::

    sage: x = PolynomialRing(RationalField(), 3, 'x').gens()
    sage: f = x[0] + x[1] - 2*x[1]*x[2]
    sage: f
    -2*x1*x2 + x0 + x1
    sage: f(1,2,0)
    3
    sage: f(1,2,5)
    -17

这也适用于有理函数：

.. link

::

    sage: h = f /(x[1] + x[2])
    sage: h
    (-2*x1*x2 + x0 + x1)/(x1 + x2)
    sage: h(1,2,3)
    -9/5

.. index::
   pair: polynomial; symbolic manipulation

Sage 还可以进行符号操作：

::

    sage: var('x,y,z')
    (x, y, z)
    sage: f = (x + 3*y + x^2*y)^3; f
    (x^2*y + x + 3*y)^3
    sage: f(x=1,y=2,z=3)
    729
    sage: f.expand()
    x^6*y^3 + 3*x^5*y^2 + 9*x^4*y^3 + 3*x^4*y + 18*x^3*y^2 +
    27*x^2*y^3 +
    x^3 + 9*x^2*y + 27*x*y^2 + 27*y^3
    sage: f(x = 5/z)
    (3*y + 25*y/z^2 + 5/z)^3
    sage: g = f.subs(x = 5/z); g
    (3*y + 25*y/z^2 + 5/z)^3
    sage: h = g.rational_simplify(); h
    (27*y^3*z^6 + 135*y^2*z^5 + 225*(3*y^3 + y)*z^4 + 125*(18*y^2 + 1)*z^3 +
    15625*y^3 + 9375*y^2*z + 1875*(3*y^3 + y)*z^2)/z^6

多元多项式的根
==============

在某些情况下，Sage（使用 Singular 接口）可以
（假设解形成零维代数簇）使用 Gröbner 基求解多元多项式方程。
以下是一个简单的示例：

::

    sage: R = PolynomialRing(QQ, 2, 'ab', order='lp')
    sage: a,b = R.gens()
    sage: I = (a^2-b^2-3, a-2*b)*R
    sage: B = I.groebner_basis(); B
    [a - 2*b, b^2 - 1]

所以 `b=\pm 1` 且 `a=2b`。

.. index:
   pair: polynomial; Groebner basis of ideal

.. _section-groebner:

Gröbner 基
==========

此计算在后台使用 Singular 来计算 Gröbner 基。

::

    sage: R = PolynomialRing(QQ, 4, 'abcd', order='lp')
    sage: a,b,c,d = R.gens()
    sage: I = (a+b+c+d, a*b+a*d+b*c+c*d, a*b*c+a*b*d+a*c*d+b*c*d, a*b*c*d-1)*R; I
    Ideal (a + b + c + d, a*b + a*d + b*c + c*d, a*b*c + a*b*d + a*c*d + b*c*d,
    a*b*c*d - 1) of Multivariate Polynomial Ring in a, b, c, d over Rational Field
    sage: B = I.groebner_basis(); B
    [a + b + c + d,
     b^2 + 2*b*d + d^2,
     b*c - b*d + c^2*d^4 + c*d - 2*d^2,
     b*d^4 - b + d^5 - d,
     c^3*d^2 + c^2*d^3 - c - d,
     c^2*d^6 - c^2*d^2 - d^4 + 1]

你可以使用多个环，而不必像在 Singular 中那样来回切换。例如，

::

    sage: a,b,c = QQ['a,b,c'].gens()
    sage: X,Y = GF(7)['X,Y'].gens()
    sage: I = ideal(a, b^2, b^3+c^3)
    sage: J = ideal(X^10 + Y^10)

    sage: I.minimal_associated_primes ()
    [Ideal (c, b, a) of Multivariate Polynomial Ring in a, b, c over Rational Field]

    sage: J.minimal_associated_primes ()     # slightly random output
    [Ideal (Y^4 + 3*X*Y^3 + 4*X^2*Y^2 + 4*X^3*Y + X^4) of Multivariate Polynomial
    Ring in X, Y over Finite Field of size 7,
     Ideal (Y^4 + 4*X*Y^3 + 4*X^2*Y^2 + 3*X^3*Y + X^4) of Multivariate Polynomial
    Ring in X, Y over Finite Field of size 7,
     Ideal (Y^2 + X^2) of Multivariate Polynomial Ring in X, Y over Finite Field
    of size 7]

所有实际工作均由 Singular 完成。

Sage 还包括 ``gfan``，它提供了计算 Gröbner 基的其他快速算法。
更多详情，请参阅参考手册中的 "Gröbner fans" 部分。
