********
初等数论
********

计算模幂
========

如何在 Sage 中计算模幂？

要想在 Sage 中计算 `51^{2006} \pmod{97}`，请输入

::

    sage: R = Integers(97)
    sage: a = R(51)
    sage: a^2006
    12

除了 ``R = Integers(97)``，你也可以输入
``R = IntegerModRing(97)``。另一种选择是使用 GMP 接口：

::

    sage: 51.powermod(99203843984,97)
    96

.. index:: discrete logs

离散对数
=============

要找到数 `x` 使得
`b^x\equiv a \pmod m` （
`a \pmod m` 的离散对数）可以使用 ``log`` 命令：

::

    sage: r = Integers(125)
    sage: b = r.multiplicative_generator()^3
    sage: a = b^17
    sage: a.log(b)
    17

这在有限域上也适用：

::

    sage: FF = FiniteField(16,"a")
    sage: a = FF.gen()
    sage: c = a^7
    sage: c.log(a)
    7

质数
====

如何在 Sage 中构造质数？

``Primes`` 类可以进行质数检测：

::

    sage: 2^(2^12)+1 in Primes()
    False
    sage: 11 in Primes()
    True

``next_prime`` 的使用一目了然：

::

    sage: next_prime(2005)
          2011

Pari 命令 ``primepi`` 通过 ``pari(x).primepi()`` 命令使用。
它返回 `\leq x` 的质数数量，例如：

::

    sage: pari(10).primepi()
          4

使用 ``primes_first_n`` 或 ``primes`` 可以检查到，
实有 `4` 个小于等于 `10` 的质数：

::

    sage: primes_first_n(5)
    [2, 3, 5, 7, 11]
    sage: list(primes(1, 10))
    [2, 3, 5, 7]

因子
====

如何在 Sage 中计算整数因子之和？

Sage 使用 ``divisors(n)`` 求 `n` 的因子列表，
使用 ``number_of_divisors(n)`` 求 `n` 的因子数量，
并使用 ``sigma(n,k)`` 求 `n` 的因子的 `k` 次幂之和
（因此 ``number_of_divisors(n)`` 和 ``sigma(n,0)`` 等价）。

例如：

::

    sage: divisors(28); sum(divisors(28)); 2*28
    [1, 2, 4, 7, 14, 28]
    56
    56
    sage: sigma(28,0); sigma(28,1); sigma(28,2)
    6
    56
    1050

.. index:: quadratic residues

二次剩余
========

尝试一下：

::

    sage: Q = quadratic_residues(23); Q
    [0, 1, 2, 3, 4, 6, 8, 9, 12, 13, 16, 18]
    sage: N = [x for x in range(22) if kronecker(x,23)==-1]; N
    [5, 7, 10, 11, 14, 15, 17, 19, 20, 21]

Q 是模 23 的二次剩余集，N 是非剩余集。

下面是使用 ``kronecker`` 命令（也称为“勒让德符号”）构造上述集合的另一种方法：

::

    sage: [x for x in range(22) if kronecker(x,23)==1]
    [1, 2, 3, 4, 6, 8, 9, 12, 13, 16, 18]
    sage: [x for x in range(22) if kronecker(x,23)==-1]
    [5, 7, 10, 11, 14, 15, 17, 19, 20, 21]
