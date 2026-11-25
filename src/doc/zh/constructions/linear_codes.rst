.. _chapter-codes:

************
线性码与密码
************

码
==

长度为 `n` 的线性码是 `GF(q)^n` 的有限维子空间。
Sage 能够在一定程度上计算线性纠错码。它基本上对 GAP 和 GUAVA 命令进行了一些封装。
GUAVA 2.8 未包含在 Sage 4.0 的 GAP 安装中，但可以作为可选包安装。

.. index:
   pair: codes; linear
   pair: codes; Hamming

Sage 能够计算汉明 (Hamming) 码

::

    sage: C = codes.HammingCode(GF(3), 3)
    sage: C
    [13, 10] Hamming Code over GF(3)
    sage: C.minimum_distance()
    3
    sage: C.generator_matrix()
    [1 0 0 0 0 0 0 0 0 0 1 2 0]
    [0 1 0 0 0 0 0 0 0 0 0 1 2]
    [0 0 1 0 0 0 0 0 0 0 1 0 2]
    [0 0 0 1 0 0 0 0 0 0 1 1 1]
    [0 0 0 0 1 0 0 0 0 0 1 1 2]
    [0 0 0 0 0 1 0 0 0 0 2 0 2]
    [0 0 0 0 0 0 1 0 0 0 1 2 1]
    [0 0 0 0 0 0 0 1 0 0 2 1 1]
    [0 0 0 0 0 0 0 0 1 0 2 2 0]
    [0 0 0 0 0 0 0 0 0 1 0 1 1]

.. index::
   pair: codes; Golay

四种戈莱 (Golay) 码

::

    sage: C = codes.GolayCode(GF(3))
    sage: C
    [12, 6, 6] Extended Golay code over GF(3)
    sage: C.minimum_distance()
    6
    sage: C.generator_matrix()
    [1 0 0 0 0 0 2 0 1 2 1 2]
    [0 1 0 0 0 0 1 2 2 2 1 0]
    [0 0 1 0 0 0 1 1 1 0 1 1]
    [0 0 0 1 0 0 1 1 0 2 2 2]
    [0 0 0 0 1 0 2 1 2 2 0 1]
    [0 0 0 0 0 1 0 2 1 2 2 1]

以及二进制 Reed-Muller 码、二次剩余码、准二次剩余码、“随机”线性码和由满秩矩阵生成的码（通常使用行作为基）。

.. index::
   pair: codes; check matrix
   pair: codes; generator matrix
   pair: codes; dual

对于给定的码 `C`，Sage 可以返回生成矩阵、检验矩阵和对偶码：

::

    sage: C = codes.HammingCode(GF(2), 3)
    sage: Cperp = C.dual_code()
    sage: C; Cperp
    [7, 4] Hamming Code over GF(2)
    [7, 3] linear code over GF(2)
    sage: C.generator_matrix()
      [1 0 0 0 0 1 1]
      [0 1 0 0 1 0 1]
      [0 0 1 0 1 1 0]
      [0 0 0 1 1 1 1]
    sage: C.parity_check_matrix()
      [1 0 1 0 1 0 1]
      [0 1 1 0 0 1 1]
      [0 0 0 1 1 1 1]
    sage: C.dual_code()
    [7, 3] linear code over GF(2)
    sage: C = codes.HammingCode(GF(4,'a'), 3)
    sage: C.dual_code()
    [21, 3] linear code over GF(4)

对于 `C` 和向量 `v\in GF(q)^n`，
Sage 可以尝试使用校验子解码 (syndrome decoding) 来解码 `v`
（即，在汉明度量下找到最接近 `v` 的码字 `c\in C`）。
到目前为止，尚未实现任何特殊的解码方法。

::

    sage: C = codes.HammingCode(GF(2), 3)
    sage: MS = MatrixSpace(GF(2),1,7)
    sage: F = GF(2); a = F.gen()
    sage: v = vector([a,a,F(0),a,a,F(0),a])
    sage: c = C.decode_to_code(v, "Syndrome"); c
    (1, 1, 0, 1, 0, 0, 1)
    sage: c in C
    True

要绘制码的权重分布（直方图），可以使用 Sage 附带的 matplotlib 包：

::

    sage: C = codes.HammingCode(GF(2), 4)
    sage: C
    [15, 11] Hamming Code over GF(2)
    sage: w = C.weight_distribution(); w
     [1, 0, 0, 35, 105, 168, 280, 435, 435, 280, 168, 105, 35, 0, 0, 1]
    sage: J = range(len(w))
    sage: W = IndexedSequence([ZZ(w[i]) for i in J],J)
    sage: P = W.plot_histogram()

现在输入 ``show(P)`` 来展示此内容。

我们完全跳过了几个编码理论函数。请参阅参考手册或文件
``coding/linear_codes.py`` 以获取示例。

Sage 还可以通过 Singular 接口 § sec:agcodes 计算代数几何码，称为 AG 码。
还可以通过 Sage 接口将 GUAVA 的 AG 码实现用于 GAP ``gap_console()``。
更多细节请参阅 GUAVA 手册。 {GUAVA}

密码
====

线性反馈移位寄存器 (LFSR)
----------------------------

在 Sage 中实现了一种特殊类型的流密码，即定义在有限域上的线性反馈移位寄存器（LFSR）序列。
流密码长期以来一直被用作伪随机数生成器的来源。
{linear feedback shift register}

S. Golomb {G} 给出了数字序列 `{\bf a}=\{a_n\}_{n=1}^\infty` (`a_n\in \{0,1\}`)
被认为是“随机”的三个统计特征。
将 `{\bf a}` 的自相关定义为：

.. MATH::

   C(k)=C(k,{\bf a})=\lim_{N\rightarrow \infty}
   \frac{1}{N}\sum_{n=1}^N (-1)^{a_n+a_{n+k}}.

如果 `a` 是周期性的，周期为 `P`，则该公式可以简化为：

.. MATH::

   C(k)=\frac{1}{P}\sum_{n=1}^P (-1)^{a_n+a_{n+k}}

假设 `a` 是周期性的，周期为 `P`。


-  平衡： `|\sum_{n=1}^P(-1)^{a_n}|\leq 1`.

-  低自相关：

   .. MATH::

      C(k)=
      \left\{
      \begin{array}{cc}
      1,& k=0,\\
      \epsilon, & k\not= 0.
      \end{array}
      \right.

   （对于满足头两个特性的序列，已知 `\epsilon=-1/P` 必然成立。）

-  比例运行特性：在每个周期中，一半的递推序列长度为 `1`, 四分之一长度为 `2`, 等等。
   此外，`1` 的递推次数和 `0` 的递推次数相等。


满足这些特性的序列将被称为伪随机数列。{pseudo-random}

一般反馈移位寄存器是一个映射 `f:{\bf F}_q^d\rightarrow {\bf F}_q^d`，形式为：

.. MATH::

   \begin{array}{c}
   f(x_0,...,x_{n-1})=(x_1,x_2,...,x_n),\\
   x_n=C(x_0,...,x_{n-1}),
   \end{array}


其中 `C:{\bf F}_q^d\rightarrow {\bf F}_q` 是一个给定的函数。当 `C` 形式为：

.. MATH::

    C(x_0,...,x_{n-1}) = c_0 x_0 + ... + c_{n-1} x_{n-1},

对于某些给定常数 `c_i\in {\bf F}_q`, 该映射被称为线性反馈移位寄存器（LFSR）。
系数序列 `c_i` 被称为密钥，多项式

.. MATH::

   C(x) = 1+ c_0x + \dots + c_{n-1}x^n

.. index::
   pair: ciphers; connection polynomial

有时被称为连接多项式。


例如：在 `GF(2)` 上，如果
`[c_0,c_1,c_2,c_3]=[1,0,0,1]` 则
`C(x) = 1 + x + x^4`,

.. MATH::

x_n = x_{n-4} + x_{n-1},\ \ \ n\geq 4.


那么 LFSR 序列为

.. MATH::

   \begin{array}{c}
   1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, \\
   1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, ...\ .
   \end{array}


此 `0,1` 序列是周期性的，周期为 `P=2^4-1=15`，并满足 Golomb 的三个随机性条件。
然而，这个以 15 为周期的序列可以通过只知道 8 个项被“破解”（即生成 `g(x)` 的过程）！
这是 Berlekamp-Massey 算法 {M} 的功能，已通过 ``lfsr_connection_polynomial`` 实现
（产生 ``berlekamp_massey`` 的反向结果）。

::

    sage: F = GF(2)
    sage: o = F(0)
    sage: l = F(1)
    sage: key = [l,o,o,l]; fill = [l,l,o,l]; n = 20
    sage: s = lfsr_sequence(key,fill,n); s
    [1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0]
    sage: lfsr_autocorrelation(s,15,7)
    4/15
    sage: lfsr_autocorrelation(s,15,0)
    8/15
    sage: lfsr_connection_polynomial(s)
    x^4 + x + 1
    sage: from sage.matrix.berlekamp_massey import berlekamp_massey
    sage: berlekamp_massey(s)
    x^4 + x^3 + 1

经典密码
--------

有一个用于加密系统的类型（由 David Kohel 创建，他还编写了下面的示例），实现了经典加密系统。一般接口如下：

::

    sage: S = AlphabeticStrings()
    sage: S
    Free alphabetic string monoid on A-Z
    sage: E = SubstitutionCryptosystem(S)
    sage: E
    Substitution cryptosystem on Free alphabetic string monoid on A-Z
    sage: K = S([ 25-i for i in range(26) ])
    sage: e = E(K)
    sage: m = S("THECATINTHEHAT")
    sage: e(m)
    GSVXZGRMGSVSZG

下面是另一个例子：

::

    sage: S = AlphabeticStrings()
    sage: E = TranspositionCryptosystem(S,15);
    sage: m = S("THECATANDTHEHAT")
    sage: G = E.key_space()
    sage: G
    Symmetric group of order 15! as a permutation group
    sage: g = G([ 3, 2, 1, 6, 5, 4, 9, 8, 7, 12, 11, 10, 15, 14, 13 ])
    sage: e = E(g)
    sage: e(m)
    EHTTACDNAEHTTAH

其思想是加密系统是一个映射 `E: KS \to \text{Hom}_\text{Set}(MS,CS)`，
其中 math:`KS`, `MS` 和 `CS` 分别是密钥空间、明文（或消息）空间和密文空间。
假设 `E` 为单射，所以 ``e.key()`` 返回原像密钥。
