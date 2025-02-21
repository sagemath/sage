.. Number Theory

数論
=============

Sageは数論関連の多彩な機能を備えている．
例えば，以下のようにして :math:`\ZZ/N\ZZ` 上の演算を実行することができる:

::

    sage: R = IntegerModRing(97)
    sage: a = R(2) / R(3)
    sage: a
    33
    sage: a.rational_reconstruction()
    2/3
    sage: b = R(47)
    sage: b^20052005
    50
    sage: b.modulus()
    97
    sage: b.is_square()
    True

Sageは数論では標準となっている関数群を装備している．例えば


::

    sage: gcd(515,2005)
    5
    sage: factor(2005)
    5 * 401
    sage: c = factorial(25); c
    15511210043330985984000000
    sage: [valuation(c,p) for p in prime_range(2,23)]
    [22, 10, 6, 3, 2, 1, 1, 1]
    sage: next_prime(2005)
    2011
    sage: previous_prime(2005)
    2003
    sage: divisors(28); sum(divisors(28)); 2*28
    [1, 2, 4, 7, 14, 28]
    56
    56

という具合で，申し分なしだ．


Sageの ``sigma(n,k)`` 関数は， ``n`` の商の :math:`k` 乗の和を計算する:

::

    sage: sigma(28,0); sigma(28,1); sigma(28,2)
    6
    56
    1050

以下では，拡張ユークリッド互除法，オイラーの :math:`\phi` -関数，そして中国剰余定理を見てみよう:


::

    sage: d,u,v = xgcd(12,15)
    sage: d == u*12 + v*15
    True
    sage: n = 2005
    sage: inverse_mod(3,n)
    1337
    sage: 3 * 1337
    4011
    sage: prime_divisors(n)
    [5, 401]
    sage: phi = n*prod([1 - 1/p for p in prime_divisors(n)]); phi
    1600
    sage: euler_phi(n)
    1600
    sage: prime_to_m_part(n, 5)
    401

次に， :math:`3n+1` 問題をちょっと調べてみる．

::

    sage: n = 2005
    sage: for i in range(1000):
    ....:     n = 3 * odd_part(n) + 1
    ....:     if odd_part(n) == 1:
    ....:         print(i)
    ....:         break
    38

最後に，中国剰余定理を確かめてみよう．


::

    sage: x = crt(2, 1, 3, 5); x
    11
    sage: x % 3  # x mod 3 = 2
    2
    sage: x % 5  # x mod 5 = 1
    1
    sage: [binomial(13,m) for m in range(14)]
    [1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1]
    sage: [binomial(13,m)%2 for m in range(14)]
    [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
    sage: [kronecker(m,13) for m in range(1,13)]
    [1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, 1]
    sage: n = 10000; sum([moebius(m) for m in range(1,n)])
    -23
    sage: Partitions(4).list()
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]


:math:`p` \-進数
------------------------

Sageには :math:`p` \-進数体も組込まれている．
ただし，いったん生成された :math:`p` \-進数体については，後でその精度を変更することはできないことを注意しておく．


::

    sage: K = Qp(11); K
    11-adic Field with capped relative precision 20
    sage: a = K(211/17); a
    4 + 4*11 + 11^2 + 7*11^3 + 9*11^5 + 5*11^6 + 4*11^7 + 8*11^8 + 7*11^9
      + 9*11^10 + 3*11^11 + 10*11^12 + 11^13 + 5*11^14 + 6*11^15 + 2*11^16
      + 3*11^17 + 11^18 + 7*11^19 + O(11^20)
    sage: b = K(3211/11^2); b
    10*11^-2 + 5*11^-1 + 4 + 2*11 + O(11^18)

:math:`p` \-進数体あるいは :math:`QQ` 以外の数体上に整数環を実装するために多大の労力が投入されてきている．
興味ある読者はGoogleグループ ``sage-support`` で専門家に詳細を聞いてみてほしい．


``NumberField`` クラスには，すでに多くの関連メソッドが実装されている．


::

    sage: R.<x> = PolynomialRing(QQ)
    sage: K = NumberField(x^3 + x^2 - 2*x + 8, 'a')
    sage: K.integral_basis()
    [1, 1/2*a^2 + 1/2*a, a^2]

.. link

::

    sage: K.galois_group()
    Galois group 3T2 (S3) with order 6 of x^3 + x^2 - 2*x + 8

.. link

::

    sage: K.polynomial_quotient_ring()
    Univariate Quotient Polynomial Ring in a over Rational Field with modulus
    x^3 + x^2 - 2*x + 8
    sage: K.units()
    (-3*a^2 - 13*a - 13,)
    sage: K.discriminant()
    -503
    sage: K.class_group()
    Class group of order 1 of Number Field in a with
    defining polynomial x^3 + x^2 - 2*x + 8
    sage: K.class_number()
    1
