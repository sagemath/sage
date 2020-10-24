.. _section-poly:

多項式
===========

この節では，Sage上で多項式を生成し利用する方法について解説する．


.. _section-univariate:


1変数多項式
----------------------

多項式環を生成するには，三通りのやり方がある．

::

    sage: R = PolynomialRing(QQ, 't')
    sage: R
    Univariate Polynomial Ring in t over Rational Field

このやり方では，多項式環を生成し，画面表示時の変数名として(文字列) ``t`` を割り当てている．
ただし，これはSage内で使うために記号 ``t`` を定義しているわけではないから，( :math:`t^2+1` のようにして)  ``R`` 上の多項式を入力するときには使えないことを注意しておく．

これに代わるやり方として

.. link

::

    sage: S = QQ['t']
    sage: S == R
    True

があるが， 記号 ``t`` については最初の方法と同じ問題が残る．

第三の，とても便利な方法が

::

    sage: R.<t> = PolynomialRing(QQ)

とするか
::

    sage: R.<t> = QQ['t']

あるいはまた

::

    sage: R.<t> = QQ[]

とすることである．

最後の方法には，多項式の変数として変数 ``t`` が定義される余禄がついてくるので，以下のようにして簡単に ``R`` の元を構成することができる．
(第三の方法が，Magmaにおけるコンストラクタ記法によく似ていることに注意．Magmaと同じように，Sageでも多様なオブジェクトでコンストラクタ記法を使うことができる．)

.. link

::

    sage: poly = (t+1) * (t+2); poly
    t^2 + 3*t + 2
    sage: poly in R
    True

どの方法で多項式環を定義していても，その変数を :math:`0` 番目の生成元として取り出すことができる:

::

    sage: R = PolynomialRing(QQ, 't')
    sage: t = R.0
    sage: t in R
    True


Sageでは，複素数も多項式環と似た構成法で生成されている．
すなわち，実数体上に記号 ``i`` を生成元として構成されているのが複素数であると見なすことができるのだ．
それを示すのが:

::

    sage: CC
    Complex Field with 53 bits of precision
    sage: CC.0  # CCの0番目の生成元
    1.00000000000000*I


多項式環を生成するときは，以下のように環と生成元の両方を同時に作るか，あるいは生成元のみを作るか選ぶことができる:

::

    sage: R, t = QQ['t'].objgen()
    sage: t    = QQ['t'].gen()
    sage: R, t = objgen(QQ['t'])
    sage: t    = gen(QQ['t'])

最後に， :math:`\QQ[t]` 上の演算を試してみよう．

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


因数分解の際には，定数(倍)部がきちんと分離して記録されていることに注目．


仕事中に，例えば ``R.cyclotomic_polynomial`` 関数を頻繁に使う必要があったとしよう．
研究発表の際には，Sage本体だけではなく，円周等分多項式の実質的な計算に使われたコンポーネントを突きとめて引用するのが筋だ．
この場合、 ``R.cyclotomic_polynomial??`` と入力してソースコードを表示すると，すぐに ``f = pari.polcyclo(n)`` という行が見つかるはずだ．
これで円周等分多項式の計算にPARIが使われていることが判ったのだから，発表の際にはPARIも引用することができる．

多項式同士で割算すると，結果は(Sageが自動的に生成する)有理数体の元になる．


::

    sage: x = QQ['x'].0
    sage: f = x^3 + 1; g = x^2 - 17
    sage: h = f/g;  h
    (x^3 + 1)/(x^2 - 17)
    sage: h.parent()
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

``QQ[x]`` の有理数体上でローラン級数による級数展開を計算することができる:

::

    sage: R.<x> = LaurentSeriesRing(QQ); R
    Laurent Series Ring in x over Rational Field
    sage: 1/(1-x) + O(x^10)
    1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)


異なる変数名を割り当てて生成した1変数多項式環は，ぞれぞれ異なる環と見なされる．


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


環は変数名によって識別される．
同じ変数 ``x``  を使うと，異なる環をもう1つ作ったつもりでいても，そうはならないことに注意してほしい．


::

    sage: R = PolynomialRing(QQ, "x")
    sage: T = PolynomialRing(QQ, "x")
    sage: R == T
    True
    sage: R is T
    True
    sage: R.0 == T.0
    True

Sageでは，任意の基底環上で巾級数環およびローラン級数環を扱うことができる．
以下の例では， :math:`\GF{7}[[T]]` の元を生成し，ついでその逆数をとって :math:`\GF{7}((T))` の元を作っている．

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

巾級数環を生成するには，二重括弧を使う省略記法を用いることもできる:


::

    sage: GF(7)[['T']]
    Power Series Ring in T over Finite Field of size 7



多変数多項式
------------------------

複数個の変数を含む多項式を扱うには，まず多項式環と変数を宣言する．
::

    sage: R = PolynomialRing(GF(5),3,"z") # 3 = 変数の数
    sage: R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5


1変数の多項式を定義したときと同じように，他の方法もある:

::

    sage: GF(5)['z0, z1, z2']
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5
    sage: R.<z0,z1,z2> = GF(5)[]; R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5


さらに，変数名を1文字にしたければ，以下のような略記法を使えばよい:

::

    sage: PolynomialRing(GF(5), 3, 'xyz')
    Multivariate Polynomial Ring in x, y, z over Finite Field of size 5


ここで，ちょっと計算してみよう．

::

    sage: z = GF(5)['z0, z1, z2'].gens()
    sage: z
    (z0, z1, z2)
    sage: (z[0]+z[1]+z[2])^2
    z0^2 + 2*z0*z1 + z1^2 + 2*z0*z2 + 2*z1*z2 + z2^2


多項式環を生成するには，もっと数学寄りの記号法を使うこともできる．


::

    sage: R = GF(5)['x,y,z']
    sage: x,y,z = R.gens()
    sage: QQ['x']
    Univariate Polynomial Ring in x over Rational Field
    sage: QQ['x,y'].gens()
    (x, y)
    sage: QQ['x'].objgens()
    (Univariate Polynomial Ring in x over Rational Field, (x,))


Sageの多変数多項式は，多項式に対する分配表現(distributive representation)とPyhonのディクショナリを使って実装されている．
gcdやイデアルのグレブナー基底の計算にはSingular [Si]_ を経由している部分がある．

::

    sage: R, (x, y) = PolynomialRing(RationalField(), 2, 'xy').objgens()
    sage: f = (x^3 + 2*y^2*x)^2
    sage: g = x^2*y^2
    sage: f.gcd(g)
    x^2


次に， :math:`f` と :math:`g` から生成されるイデアル :math:`(f,g)` を求めてみる．
これには ``(f,g)`` に ``R`` を掛けてやるだけでよい(``ideal([f,g])`` あるいは ``ideal(f,g)`` としても同じだ)．

.. link


::

    sage: I = (f, g)*R; I
    Ideal (x^6 + 4*x^4*y^2 + 4*x^2*y^4, x^2*y^2) of Multivariate Polynomial
    Ring in x, y over Rational Field
    sage: B = I.groebner_basis(); B
    [x^6, x^2*y^2]
    sage: x^2 in I
    False


ちなみに，上のグレブナー基底はリストではなく不変性シーケンスとして与えられている．
これは，基底がユニバースまたは親クラスとなっているため変更できないことを意味している(グレブナー基底が変えられてしまうとその基底系に依存するルーチン群も働かなくなるから当然のことだ)．

.. link


::

    sage: B.parent()
    <class 'sage.rings.polynomial.multi_polynomial_sequence.PolynomialSequence_generic'>
    sage: B.universe()
    Multivariate Polynomial Ring in x, y over Rational Field
    sage: B[1] = x
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

(種類はまだ十分ではないものの)可換代数の中にはSigular経由でSage上に実装されているものがある．
例えば， :math:`I`: の準素分解および随伴素イデアルを求めることができる:


.. link

::

    sage: I.primary_decomposition()
    [Ideal (x^2) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y^2, x^6) of Multivariate Polynomial Ring in x, y over Rational Field]
    sage: I.associated_primes()
    [Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field]
