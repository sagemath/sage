.. -*- coding: utf-8 -*-

.. _section-coercion:

================================
父结构、转换与强制转换
================================

这一节可能比前一节更技术化，但为了有效且高效地使用 Sage 中的环和其他代数结构，
理解父结构和强制转换的意义非常重要。

请注意，我们在这里只解释概念，不展示具体实现。
面向实现的教程可以参见
`Sage thematic tutorial <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_ 。

元素
--------

如果想在 Python 中实现一个环，第一步是创建一个类来表示该环的元素 ``X``，
并为其提供必要的双下划线方法，例如 ``__add__``, ``__sub__``, ``__mul__``，
同时确保环公理成立。

由于 Python 是一种强类型（但动态类型）语言，可能会想到为每个环实现一个 Python 类。
毕竟，Python 有整数类型 ``<int>`` 和实数类型 ``<float>`` 等等。
但这种方法很快就会失败：环的数量是无限的，无法实现无限多个类。

相反，可以创建一个类层次结构来实现常见的代数结构元素，例如群、环、斜域、交换环、域、代数等等。

但这意味着不同环的元素可以具有相同的类型。

::

    sage: P.<x,y> = GF(3)[]
    sage: Q.<a,b> = GF(4,'z')[]
    sage: type(x)==type(a)
    True

另一方面，也可以有不同的 Python 类来实现相同的数学结构（例如稠密矩阵与稀疏矩阵）

::

    sage: P.<a> = PolynomialRing(ZZ)
    sage: Q.<b> = PolynomialRing(ZZ, sparse=True)
    sage: R.<c> = PolynomialRing(ZZ, implementation='NTL')
    sage: type(a); type(b); type(c)
    <class 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
    <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_integral_domain_with_category.element_class'>
    <class 'sage.rings.polynomial.polynomial_integer_dense_ntl.Polynomial_integer_dense_ntl'>

这带来了两个问题：一方面，如果两个元素是相同类的实例，可以预期它们的 ``__add__`` 方法能够相加；
但如果这些元素属于非常不同的环，则不希望如此。另一方面，如果两个元素属于同一环的不同实现，想要相加，
但如果它们属于不同的 Python 类，这并不容易实现。

解决这些问题的方法称为“强制转换”，将在下面解释。

然而，每个元素都必须知道它属于哪个父结构。这可以通过 ``parent()`` 方法获得：

.. link

::

    sage: a.parent(); b.parent(); c.parent()
    Univariate Polynomial Ring in a over Integer Ring
    Sparse Univariate Polynomial Ring in b over Integer Ring
    Univariate Polynomial Ring in c over Integer Ring (using NTL)


父结构与范畴
----------------------

与代数结构元素的 Python 类层次结构类似，Sage 也提供包含这些元素的代数结构的类。
在 Sage 中包含元素的结构称为“父结构”，并且有一个基类。
大致上与数学概念的层次结构一致，有一系列类，例如集合、环、域等等：

::

    sage: isinstance(QQ,Field)
    True
    sage: isinstance(QQ, Ring)
    True
    sage: isinstance(ZZ,Field)
    False
    sage: isinstance(ZZ, Ring)
    True

在代数中，共享相同代数结构的对象被归类到所谓的“范畴”中。
因此，Sage 中类层次结构与范畴层次结构之间有一个粗略的类比。
然而，不应过分强调 Python 类与范畴的类比。毕竟，数学范畴也在 Sage 中实现：

::

    sage: Rings()
    Category of rings
    sage: ZZ.category()
    Join of Category of Dedekind domains
        and Category of euclidean domains
        and Category of noetherian rings
        and Category of infinite enumerated sets
        and Category of metric spaces
    sage: ZZ.category().is_subcategory(Rings())
    True
    sage: ZZ in Rings()
    True
    sage: ZZ in Fields()
    False
    sage: QQ in Fields()
    True

虽然 Sage 的类层次结构集中在实现细节上，但 Sage 的范畴框架更集中在数学结构上。
可以在范畴中实现不依赖具体实现的通用方法和测试。

Sage 中的父结构应该是唯一的 Python 对象。
例如，一旦创建了一个具有特定基环和特定生成器列表的多项式环，结果将被缓存：

::

    sage: RR['x','y'] is RR['x','y']
    True


类型与父结构
--------------------
类型 ``RingElement`` 并不完全对应于数学概念中的环元素。
例如，虽然方阵属于一个环，但它们不是 ``RingElement`` 的实例：

::

    sage: M = Matrix(ZZ,2,2); M
    [0 0]
    [0 0]
    sage: isinstance(M, RingElement)
    False

虽然在 Sage 中 *父结构* 是唯一的，但在一个父结构中的相等元素不一定是相同的。
这与 Python 对某些（虽然不是全部）整数的行为形成对比：

::

    sage: int(1) is int(1) # Python int
    True
    sage: int(-15) is int(-15)
    False
    sage: 1 is 1           # Sage Integer
    False

不同环的元素通常不是通过它们的类型区分，而是通过它们的父结构区分：

::

    sage: a = GF(2)(1)
    sage: b = GF(5)(1)
    sage: type(a) is type(b)
    True
    sage: parent(a)
    Finite Field of size 2
    sage: parent(b)
    Finite Field of size 5

因此，从代数的角度来看，**元素的父结构比它的类型更重要。**

转换与强制转换
--------------------------

在某些情况下，可以将一个父结构的元素转换为另一个父结构的元素。
这样的转换可以是显式的也可以是隐式的（被称为 *强制转换*）。

读者可能知道例如 C 语言中的 *类型转换* 和 *类型强制转换* 的概念。
Sage 中也有转换和强制转换的概念。但 Sage 中的概念集中在 *父结构* 上，而不是类型上。
所以请不要将 C 语言中的类型转换与 Sage 中的转换混淆！

我们在这里给出一个相当简短的说明。
详细描述和实现信息，请参阅参考手册中的强制转换章节以及
`thematic tutorial <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_.

关于在 *不同* 环的元素上进行算术运算的可能性，有两种极端观点：

* 不同的环是不同的世界，对不同环的元素进行加法或乘法没有任何意义；
  即使 ``1 + 1/2`` 也没有意义，因为第一个加数是整数，第二个是有理数。

或者

* 如果一个环 ``R1`` 的元素 ``r1`` 可以以某种方式在另一个环 ``R2`` 中解释，
  那么所有涉及 ``r1`` 和任意 ``R2`` 元素的算术运算都是允许的。
  乘法单位存在于所有域和许多环，它们应该都是相等的。

Sage 选择了一种折衷方案。如果 ``P1`` 和 ``P2`` 是父结构，``p1`` 是 ``P1`` 的元素，
那么用户可以显式请求将 ``p1`` 在 ``P2`` 中解释。这在所有情况下可能没有意义，
或者对于 ``P1`` 的所有元素都没有定义，用户需要确保其合理性。我们称之为 **转换**：

::

    sage: a = GF(2)(1)
    sage: b = GF(5)(1)
    sage: GF(5)(a) == b
    True
    sage: GF(2)(b) == a
    True

然而，只有当这种转换可以彻底和一致地完成时，才会发生 *隐式* （或自动）转换。
数学的严谨性在这一点上至关重要。

这种隐式转换称为 **强制转换**。如果定义了强制转换，那么它必须与转换一致。
定义强制转换需要满足两个条件：

#. 从 ``P1`` 到 ``P2`` 的强制转换必须由结构保持映射给出（例如环同态）。
   仅仅一些 ``P1`` 的元素可以映射到 ``P2`` 是不够的，映射必须尊重 ``P1`` 的代数结构。
#. 这些强制转换映射的选择必须一致：如果 ``P3`` 是第三个父结构，
   那么从 ``P1`` 到 ``P2`` 的选定强制转换与从 ``P2`` 到 ``P3`` 的强制转换的组合
   必须与从 ``P1`` 到 ``P3`` 的选定强制转换一致。特别是，
   如果存在从 ``P1`` 到 ``P2`` 和从 ``P2`` 到 ``P1`` 的强制转换，则组合必须是 ``P1`` 的恒等映射。

因此，尽管可以将 ``GF(2)`` 的每个元素转换为 ``GF(5)``，但不能强制转换，
因为 ``GF(2)`` 和 ``GF(5)`` 之间没有环同态。

一致性方面更难解释。我们用多元多项式环来说明。在应用中，保留名称的强制转换最有意义。因此，我们有：

::

    sage: R1.<x,y> = ZZ[]
    sage: R2 = ZZ['y','x']
    sage: R2.has_coerce_map_from(R1)
    True
    sage: R2(x)
    x
    sage: R2(y)
    y
    sage: R2.coerce(y)
    y

如果没有保留名称的环同态，则不定义强制转换。然而，转换可能仍然是可能的，即通过根据生成器列表中的位置映射环生成器：

.. link

::

    sage: R3 = ZZ['z','x']
    sage: R3.has_coerce_map_from(R1)
    False
    sage: R3(x)
    z
    sage: R3(y)
    x
    sage: R3.coerce(y)
    Traceback (most recent call last):
    ...
    TypeError: no canonical coercion
    from Multivariate Polynomial Ring in x, y over Integer Ring
    to Multivariate Polynomial Ring in z, x over Integer Ring

但这种保留位置的转换不符合强制转换：通过组合从 ``ZZ['x','y']`` 到 ``ZZ['y','x']`` 的保留名称映射
与从 ``ZZ['y','x']`` 到 ``ZZ['a','b']`` 的保留位置映射，将得到一个既不保留名称也不保留位置的映射，违反了一致性。

如果存在强制转换，它将用于比较不同环的元素或进行算术运算。这通常很方便，
但用户应该意识将 ``==`` 关系扩展到不同父结构的边界可能很容易导致过度使用。
例如，虽然 ``==`` 应该是 **同一** 环元素上的等价关系，但如果涉及 *不同* 环，则不一定如此。
例如，``ZZ`` 和有限域中的 ``1`` 被认为是相等的，因为从整数到任何有限域都有一个规范的强制转换。
然而，通常两个不同的有限域之间没有强制转换。因此我们有：

.. link

::

    sage: GF(5)(1) == 1
    True
    sage: 1 == GF(2)(1)
    True
    sage: GF(5)(1) == GF(2)(1)
    False
    sage: GF(5)(1) != GF(2)(1)
    True

同理，我们有：

.. link

::

    sage: R3(R1.1) == R3.1
    True
    sage: R1.1 == R3.1
    False
    sage: R1.1 != R3.1
    True


一致性条件的另一个结果是强制转换只能从精确环（例如有理数 ``QQ``）到不精确环（例如具有固定精度的实数 ``RR``），而不能反过来。
原因是从 ``QQ`` 到 ``RR`` 的强制转换与从 ``RR`` 到 ``QQ`` 的转换的组合应该是 ``QQ`` 上的恒等映射。
但这是不可能的，因为在 ``RR`` 中一些不同的有理数可能被视为相等，如下例所示：

::

    sage: RR(1/10^200+1/10^100) == RR(1/10^100)
    True
    sage: 1/10^200+1/10^100 == 1/10^100
    False


当比较两个父结构 ``P1`` 和 ``P2`` 的元素时，可能没有两个环之间的强制转换，
但有一个规范的父结构 ``P3`` 可选，使得 ``P1`` 和 ``P2`` 都强制转换到 ``P3``。
在这种情况下，也会发生强制转换。一个典型用例是有理数和具有整数系数的多项式之和，产生具有有理系数的多项式：

::

    sage: P1.<x> = ZZ[]
    sage: p = 2*x+3
    sage: q = 1/2
    sage: parent(p)
    Univariate Polynomial Ring in x over Integer Ring
    sage: parent(p+q)
    Univariate Polynomial Ring in x over Rational Field

注意，原则上结果在 ``ZZ['x']`` 的分数域中也有意义。
然而，Sage 会尝试选择一个 *规范的* 共同父结构，使得看起来最自然（在我们的例子中是 ``QQ['x']``）。
如果几个潜在的共同父结构看起来同样自然，为了获得可靠的结果，Sage *不会* 随机选择其中一个。
该选择所基于的机制在
`thematic tutorial <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_
中进行了解释。

以下示例不会发生强制转换到共同父结构：

::

    sage: R.<x> = QQ[]
    sage: S.<y> = QQ[]
    sage: x+y
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +: 'Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

原因是 Sage 不会选择潜在候选结构 ``QQ['x']['y']``, ``QQ['y']['x']``, ``QQ['x','y']`` 或 ``QQ['y','x']`` 之一，
因为所有这四个成对不同的结构看起来都是自然的共同父结构，并且没有明显的规范选择。
