****
数域
****

分歧 (Ramification)
===================

如何在 Sage 中计算具有给定判别式和分歧的数域？

Sage 可以访问 Jones 数域数据库，该数据库包含有界分歧且度数不超过 6 的数域。
该数据库必须单独安装（``database_jones_numfield``）。

.. index::
   pair: number field; database

首先加载数据库：

::

    sage: J = JonesDatabase()            # optional - database
    sage: J                              # optional - database
    John Jones's table of number fields with bounded ramification and degree <= 6

.. index::
   pair: number field; discriminant

列出数据库中所有分歧最多为 2 的数域的度数和判别式：

.. link

::

    sage: [(k.degree(), k.disc()) for k in J.unramified_outside([2])] # optional - database
    [(4, -2048), (2, 8), (4, -1024), (1, 1), (4, 256), (2, -4), (4, 2048), (4, 512), (4, 2048), (2, -8), (4, 2048)]

列出在 2 之外无分歧且度数恰好为 2 的域的判别式：

.. link

::

    sage: [k.disc() for k in J.unramified_outside([2],2)] # optional - database
    [8, -4, -8]

列出数据库中在 3 和 5 处有分歧的立方域的判别式：

.. link

::

    sage: [k.disc() for k in J.ramified_at([3,5],3)] # optional - database
    [-6075, -6075, -675, -135]
    sage: factor(6075)
    3^5 * 5^2
    sage: factor(675)
    3^3 * 5^2
    sage: factor(135)
    3^3 * 5

列出所有在 101 处有分歧的域：

.. link

::

    sage: J.ramified_at(101)                     # optional - database
    [Number Field in a with defining polynomial x^2 - 101,
     Number Field in a with defining polynomial x^4 - x^3 + 13*x^2 - 19*x + 361,
     Number Field in a with defining polynomial x^5 - x^4 - 40*x^3 - 93*x^2 - 21*x + 17,
     Number Field in a with defining polynomial x^5 + x^4 - 6*x^3 - x^2 + 18*x + 4,
     Number Field in a with defining polynomial x^5 + 2*x^4 + 7*x^3 + 4*x^2 + 11*x - 6]

.. index::
   pair: number field; class_number

类数 (Class numbers)
====================

如何在 Sage 中计算数域的类数 (class number)？

``class_number`` 是一个与 QuadraticField 对象相关的方法：

::

    sage: K = QuadraticField(29, 'x')
    sage: K.class_number()
    1
    sage: K = QuadraticField(65, 'x')
    sage: K.class_number()
    2
    sage: K = QuadraticField(-11, 'x')
    sage: K.class_number()
    1
    sage: K = QuadraticField(-15, 'x')
    sage: K.class_number()
    2
    sage: K.class_group()
    Class group of order 2 with structure C2 of Number Field in x with defining polynomial x^2 + 15 with x = 3.872983346207417?*I
    sage: K = QuadraticField(401, 'x')
    sage: K.class_group()
    Class group of order 5 with structure C5 of Number Field in x with defining polynomial x^2 - 401 with x = 20.02498439450079?
    sage: K.class_number()
    5
    sage: K.discriminant()
    401
    sage: K = QuadraticField(-479, 'x')
    sage: K.class_group()
    Class group of order 25 with structure C25 of Number Field in x with defining polynomial x^2 + 479 with x = 21.88606862823929?*I
    sage: K.class_number()
    25
    sage: K.pari_polynomial()
    x^2 + 479
    sage: K.degree()
    2

下面是一个更为一般的数域类型的例子：

::

    sage: x = PolynomialRing(QQ, 'x').gen()
    sage: K = NumberField(x^5+10*x+1, 'a')
    sage: K
    Number Field in a with defining polynomial x^5 + 10*x + 1
    sage: K.degree()
    5
    sage: K.pari_polynomial()
    x^5 + 10*x + 1
    sage: K.discriminant()
    25603125
    sage: K.class_group()
    Class group of order 1 of Number Field in a with defining
    polynomial x^5 + 10*x + 1
    sage: K.class_number()
    1


-  另请参见 Math World 网站上的类数链接
   http://mathworld.wolfram.com/ClassNumber.html
   获取表格、公式和背景信息。

.. index::
   pair: number field; cyclotomic

-  对于循环域，可以尝试：

   ::

       sage: K = CyclotomicField(19)
       sage: K.class_number()    # long time
       1


更多详情，请参见 ``ring/number_field.py`` 中的文档。

.. index::
   pair: number field; integral basis

整基 (Integral basis)
=====================

如何在 Sage 中计算数域的整基？

Sage 可以计算数域的元素列表，该列表是该数域的整数全环的基。

::

    sage: x = PolynomialRing(QQ, 'x').gen()
    sage: K = NumberField(x^5+10*x+1, 'a')
    sage: K.integral_basis()
    [1, a, a^2, a^3, a^4]

接下来我们计算立方域的整数环，其中 2 是“基本判别式因子”，因此整数环不是由单个元素生成的。

::

    sage: x = PolynomialRing(QQ, 'x').gen()
    sage: K = NumberField(x^3 + x^2 - 2*x + 8, 'a')
    sage: K.integral_basis()
    [1, 1/2*a^2 + 1/2*a, a^2]
