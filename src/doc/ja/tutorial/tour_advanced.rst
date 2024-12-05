
より進んだ数学
==============================


代数幾何
------------------

Sageでは，任意の代数多様体を定義することができるが，その非自明な機能は :math:`\QQ` 上の環あるいは有限体でしか使えない場合がある．
例として，2本のアフィン平面曲線の和を取り，ついで元の曲線を和の既約成分として分離してみよう．


::

    sage: x, y = AffineSpace(2, QQ, 'xy').gens()
    sage: C2 = Curve(x^2 + y^2 - 1)
    sage: C3 = Curve(x^3 + y^3 - 1)
    sage: D = C2 + C3
    sage: D
    Affine Plane Curve over Rational Field defined by
       x^5 + x^3*y^2 + x^2*y^3 + y^5 - x^3 - y^3 - x^2 - y^2 + 1
    sage: D.irreducible_components()
    [
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x^2 + y^2 - 1,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x^3 + y^3 - 1
    ]

以上の2本の曲線の交わりを取れば，全ての交点を求めてその既約成分を計算することもできる．


.. link

::

    sage: V = C2.intersection(C3)
    sage: V.irreducible_components()
    [
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      y - 1,
      x,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      y,
      x - 1,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x + y + 2,
      2*y^2 + 4*y + 3
    ]

というわけで，点 :math:`(1,0)` および :math:`(0,1)` が双方の曲線上にあるのはすぐ見てとることができるし，
:math:`y` 成分が :math:`2y^2 + 4y + 3=0` を満足する(2次の)点についても同じことだ．


Sageでは，3次元射影空間における捻れ3次曲線のトーリック・イデアルを計算することができる:


::

    sage: R.<a,b,c,d> = PolynomialRing(QQ, 4)
    sage: I = ideal(b^2-a*c, c^2-b*d, a*d-b*c)
    sage: F = I.groebner_fan(); F
    Groebner fan of the ideal:
    Ideal (b^2 - a*c, c^2 - b*d, -b*c + a*d) of Multivariate Polynomial Ring
    in a, b, c, d over Rational Field
    sage: F.reduced_groebner_bases ()
    [[-c^2 + b*d, -b*c + a*d, -b^2 + a*c],
     [-b*c + a*d, -c^2 + b*d, b^2 - a*c],
     [-c^3 + a*d^2, -c^2 + b*d, b*c - a*d, b^2 - a*c],
     [-c^2 + b*d, b^2 - a*c, b*c - a*d, c^3 - a*d^2],
     [-b*c + a*d, -b^2 + a*c, c^2 - b*d],
     [-b^3 + a^2*d, -b^2 + a*c, c^2 - b*d, b*c - a*d],
     [-b^2 + a*c, c^2 - b*d, b*c - a*d, b^3 - a^2*d],
     [c^2 - b*d, b*c - a*d, b^2 - a*c]]
    sage: F.polyhedralfan()
    Polyhedral fan in 4 dimensions of dimension 4


楕円曲線
---------------

Sageの楕円曲線部門にはPARIの楕円曲線機能の大部分が取り込まれており，Cremonaの管理するオンラインデータベースに接続することもできる(これにはデータベースパッケージを追加する必要がある)．
さらに、Second-descentによって楕円曲線の完全Mordell-Weil群を計算するmwrankの機能が使えるし，SEAアルゴリズムの実行や同種写像全ての計算なども可能だ．
:math:`\QQ` 上の曲線群を扱うためのコードは大幅に更新され，Denis Simonによる代数的降下法ソフトウェアも取り込まれている．


楕円曲線を生成するコマンド ``EllipticCurve`` には，さまざまな書法がある:


-  EllipticCurve([:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`, :math:`a_6` ]):
   楕円曲線

   .. math::  y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6,

   を生成する．
   ただし :math:`a_i` は :math:`a_1` のペアレントクラスに合わせて型強制される．
   全ての :math:`a_i` がペアレント :math:`\ZZ` を持つ場合， :math:`a_i` は :math:`\QQ` に型強制される．


-  EllipticCurve([:math:`a_4`, :math:`a_6` ]): :math:`a_1=a_2=a_3=0` となる以外は上と同じ．


-  EllipticCurve(ラベル): Cremonaの(新しい)分類ラベルを指定して，Cremonaデータベースに登録された楕円曲線を生成する．
   ラベルは    ``"11a"`` や ``"37b2"`` といった文字列で，(以前のラベルと混同しないように)小文字でなければならない．


-  EllipticCurve(j): :math:`j` -不変量 :math:`j` を持つ楕円曲線を生成する．


-  EllipticCurve(R,[:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`,  :math:`a_6` ]):
   最初と同じように :math:`a_i` を指定して環 :math:`R` 上の楕円曲線を生成する．


以上の各コンストラクタを実際に動かしてみよう:


::

    sage: EllipticCurve([0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve([GF(5)(0),0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

    sage: EllipticCurve([1,2])
    Elliptic Curve defined by y^2  = x^3 + x + 2 over Rational Field

    sage: EllipticCurve('37a')
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve_from_j(1)
    Elliptic Curve defined by y^2 + x*y = x^3 + 36*x + 3455 over Rational Field

    sage: EllipticCurve(GF(5), [0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

点 :math:`(0,0)` は、 :math:`y^2 + y = x^3 - x` で定義される楕円曲線 :math:`E` 上にある．
Sageを使ってこの点を生成するには， ``E([0,0])`` と入力する．
Sageは，そうした楕円曲線上に点を付け加えていくことができる(楕円曲線は，無限遠点が零元、同一曲線上の3点を加えると0となる加法群としての構造を備えている):

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
    sage: P = E([0,0])
    sage: P + P
    (1 : 0 : 1)
    sage: 10*P
    (161/16 : -2065/64 : 1)
    sage: 20*P
    (683916417/264517696 : -18784454671297/4302115807744 : 1)
    sage: E.conductor()
    37

複素数体上の楕円曲線は， :math:`j` -不変量によって記述される．
Sageでは， :math:`j` -不変量を以下のようにして計算する:

::

    sage: E = EllipticCurve([0,0,0,-4,2]); E
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: E.conductor()
    2368
    sage: E.j_invariant()
    110592/37

:math:`E` と同じ :math:`j` -不変量を指定して楕円曲線を作っても，それが :math:`E` と同型になるとは限らない．
次の例でも，２つの曲線は導手(conductor)が異なるため同型にならない．


::

    sage: F = EllipticCurve_from_j(110592/37)
    sage: F.conductor()
    37

しかし， :math:`F` を2で捻ったツイスト(twist)は同型の曲線になる．


.. link

::

    sage: G = F.quadratic_twist(2); G
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: G.conductor()
    2368
    sage: G.j_invariant()
    110592/37

楕円曲線に随伴する :math:`L` -級数，あるいはモジュラー形式 :math:`\sum_{n=0}^\infty a_nq^n` の係数 :math:`a_n` を求めることもできる．
計算にはPARIのC-ライブラリを援用している:

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E.anlist(30)
    [0, 1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0, -4,
     3, 10, 2, 0, -1, 4, -9, -2, 6, -12]
    sage: v = E.anlist(10000)

:math:`a_n` を :math:`n\leq 10^5` の全てについて計算しても1秒ほどしかかからない:


.. skip

::

    sage: %time v = E.anlist(100000)
    CPU times: user 0.98 s, sys: 0.06 s, total: 1.04 s
    Wall time: 1.06


楕円曲線を，対応するCremonaの分類ラベルを指定して生成する方法もある．
そうすると，目的の楕円曲線がその階数，玉河数，単数基準(regulator)などの情報と共にプレロードされる:


::

    sage: E = EllipticCurve("37b2")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational
    Field
    sage: E = EllipticCurve("389a")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x  over Rational Field
    sage: E.rank()
    2
    sage: E = EllipticCurve("5077a")
    sage: E.rank()
    3

Cremonaのデータベースへ直接にアクセスすることも可能だ．


::

    sage: db = sage.databases.cremona.CremonaDatabase()
    sage: db.curves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1], 'b1': [[0, 1, 1, -23, -50], 0, 3]}
    sage: db.allcurves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1],
     'b1': [[0, 1, 1, -23, -50], 0, 3],
     'b2': [[0, 1, 1, -1873, -31833], 0, 1],
     'b3': [[0, 1, 1, -3, 1], 0, 3]}


この方法でデータベースから引き出されるデータは，むろん ``EllipticCurve`` 型のオブジェクトにはならない．
複数のフィールドから構成されたデータベースのレコードであるにすぎない．
デフォルトでSageに付属しているのは，導手が :math:`\leq 10000` の楕円曲線の情報要約からなる，Cremonaのデータベースの小型版である．
オプションで大型版のデータベースも用意されていて，こちらは導手が :math:`120000` までの全ての楕円曲線群の詳細情報を含む(2005年10月時点)．
さらに、Sage用の大規模版データベースパッケージ(2GB)では，Stein-Watkinsデータベース上の数千万種の楕円曲線を利用することができる．


ディリクレ指標
--------------------

ディリクレ指標とは，
環 :math:`R` に対する準同型写像 :math:`(\ZZ/N\ZZ)^* \to R^*` を， :math:`\gcd(N,x)>1` なる整数 :math:`x` を0と置くことによって写像
:math:`\ZZ \to R` へ拡張したものである．


::

    sage: G = DirichletGroup(12)
    sage: G.list()
    [Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1,
    Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1]
    sage: G.gens()
    (Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1)
    sage: len(G)
    4

ディリクレ群を作成したので、次にその元を一つ取って演算に使ってみよう．

.. link

::

    sage: G = DirichletGroup(21)
    sage: chi = G.1; chi
    Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6
    sage: chi.values()
    [0, 1, zeta6 - 1, 0, -zeta6, -zeta6 + 1, 0, 0, 1, 0, zeta6, -zeta6, 0, -1,
     0, 0, zeta6 - 1, zeta6, 0, -zeta6 + 1, -1]
    sage: chi.conductor()
    7
    sage: chi.modulus()
    21
    sage: chi.order()
    6
    sage: chi(19)
    -zeta6 + 1
    sage: chi(40)
    -zeta6 + 1

この指標に対してガロワ群 :math:`\text{Gal}(\QQ(\zeta_N)/\QQ)` がどう振る舞うか計算したり，法(modulus)の因数分解に相当する直積分解を実行することも可能だ．

.. link

::

    sage: chi.galois_orbit()
    [Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> -zeta6 + 1,
    Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6]

    sage: go = G.galois_orbits()
    sage: [len(orbit) for orbit in go]
    [1, 2, 2, 1, 1, 2, 2, 1]

    sage: G.decomposition()
    [
    Group of Dirichlet characters modulo 3 with values in Cyclotomic Field of order 6 and degree 2,
    Group of Dirichlet characters modulo 7 with values in Cyclotomic Field of order 6 and degree 2
    ]

次に，mod 20，ただし値が :math:`\QQ(i)` 上に収まるディリクレ指標の群を作成する:

::

    sage: K.<i> = NumberField(x^2+1)
    sage: G = DirichletGroup(20,K)
    sage: G
    Group of Dirichlet characters modulo 20 with values in Number Field in i with defining polynomial x^2 + 1


ついで， ``G`` の不変量をいくつか計算してみよう:

.. link

::

    sage: G.gens()
    (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1,
     Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> i)

    sage: G.unit_gens()
    (11, 17)
    sage: G.zeta()
    i
    sage: G.zeta_order()
    4

以下の例では、数体上でディリクレ指標を生成する．1の累乗根については、 ``DirichletGroup`` の3番目の引数として明示的に指定している．

::

    sage: x = polygen(QQ, 'x')
    sage: K = NumberField(x^4 + 1, 'a'); a = K.0
    sage: b = K.gen(); a == b
    True
    sage: K
    Number Field in a with defining polynomial x^4 + 1
    sage: G = DirichletGroup(5, K, a); G
    Group of Dirichlet characters modulo 5 with values in the group of order 8 generated by a in Number Field in a with defining polynomial x^4 + 1
    sage: chi = G.0; chi
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> a^2
    sage: [(chi^i)(2) for i in range(4)]
    [1, a^2, -1, -a^2]


ここで ``NumberField(x^4 + 1, 'a')`` と指定したのは，Sageに記号 `a` を使って ``K`` の内容(`a` で生成される数体上の多項式 :math:`x^4 + 1`)を表示させるためである．
その時点で記号名 `a` はいったん未定義になるが、 ``a = K.0`` (``a = K.gen()`` としても同じ)が実行されると記号 `a` は多項式 :math:`x^4+1` の根を表すようになる．


モジュラー形式
-----------------

Sageを使ってモジュラー空間の次元，モジュラー・シンポルの空間，Hecke演算子、素因数分解などを含むモジュラー形式に関連した計算を実行することができる．

モジュラー形式が張る空間の次元を求める関数が数種類用意されている．
例えば


::

    sage: from sage.modular.dims import dimension_cusp_forms
    sage: dimension_cusp_forms(Gamma0(11),2)
    1
    sage: dimension_cusp_forms(Gamma0(1),12)
    1
    sage: dimension_cusp_forms(Gamma1(389),2)
    6112

次に、レベル :math:`1` ，ウェイト :math:`12` のモジュラー・シンボル空間上でHecke演算子を計算してみよう．


::

    sage: M = ModularSymbols(1,12)
    sage: M.basis()
    ([X^8*Y^2,(0,0)], [X^9*Y,(0,0)], [X^10,(0,0)])
    sage: t2 = M.T(2)
    sage: t2
    Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1)
    of weight 12 with sign 0 over Rational Field
    sage: t2.matrix()
    [ -24    0    0]
    [   0  -24    0]
    [4860    0 2049]
    sage: f = t2.charpoly('x'); f
    x^3 - 2001*x^2 - 97776*x - 1180224
    sage: factor(f)
    (x - 2049) * (x + 24)^2
    sage: M.T(11).charpoly('x').factor()
    (x - 285311670612) * (x - 534612)^2

:math:`\Gamma_0(N)` と :math:`\Gamma_1(N)` の空間を生成することもできる．


::

    sage: ModularSymbols(11,2)
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
     0 over Rational Field
    sage: ModularSymbols(Gamma1(11),2)
    Modular Symbols space of dimension 11 for Gamma_1(11) of weight 2 with
    sign 0 over Rational Field

特性多項式と :math:`q` -展開を計算してみよう．


::

    sage: M = ModularSymbols(Gamma1(11),2)
    sage: M.T(2).charpoly('x')
    x^11 - 8*x^10 + 20*x^9 + 10*x^8 - 145*x^7 + 229*x^6 + 58*x^5 - 360*x^4
         + 70*x^3 - 515*x^2 + 1804*x - 1452
    sage: M.T(2).charpoly('x').factor()
    (x - 3) * (x + 2)^2 * (x^4 - 7*x^3 + 19*x^2 - 23*x + 11)
            * (x^4 - 2*x^3 + 4*x^2 + 2*x + 11)
    sage: S = M.cuspidal_submodule()
    sage: S.T(2).matrix()
    [-2  0]
    [ 0 -2]
    sage: S.q_expansion_basis(10)
    [
        q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)
    ]

モジュラー・シンボルの空間を，指標を指定して生成することも可能だ．

::

    sage: G = DirichletGroup(13)
    sage: e = G.0^2
    sage: M = ModularSymbols(e,2); M
    Modular Symbols space of dimension 4 and level 13, weight 2, character
    [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
    sage: M.T(2).charpoly('x').factor()
    (x - zeta6 - 2) * (x - 2*zeta6 - 1) * (x + zeta6 + 1)^2
    sage: S = M.cuspidal_submodule(); S
    Modular Symbols subspace of dimension 2 of Modular Symbols space of
    dimension 4 and level 13, weight 2, character [zeta6], sign 0, over
    Cyclotomic Field of order 6 and degree 2
    sage: S.T(2).charpoly('x').factor()
    (x + zeta6 + 1)^2
    sage: S.q_expansion_basis(10)
    [
    q + (-zeta6 - 1)*q^2 + (2*zeta6 - 2)*q^3 + zeta6*q^4 + (-2*zeta6 + 1)*q^5
      + (-2*zeta6 + 4)*q^6 + (2*zeta6 - 1)*q^8 - zeta6*q^9 + O(q^10)
    ]

以下の例では，モジュラー形式によって張られる空間に対するHecke演算子の作用を，Sageでどうやって計算するかを示す．


::

    sage: T = ModularForms(Gamma0(11),2)
    sage: T
    Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of
    weight 2 over Rational Field
    sage: T.degree()
    2
    sage: T.level()
    11
    sage: T.group()
    Congruence Subgroup Gamma0(11)
    sage: T.dimension()
    2
    sage: T.cuspidal_subspace()
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
    Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: T.eisenstein_subspace()
    Eisenstein subspace of dimension 1 of Modular Forms space of dimension 2
    for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: M = ModularSymbols(11); M
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
    0 over Rational Field
    sage: M.weight()
    2
    sage: M.basis()
    ((1,0), (1,8), (1,9))
    sage: M.sign()
    0

:math:`T_p` は通常のHecke演算子( :math:`p` は素数)を表す．
Hecke演算子 :math:`T_2` ， :math:`T_3` ， :math:`T_5` はモジュラー・シンボル空間にどんな作用を及ぼすのだろうか？


.. link

::

    sage: M.T(2).matrix()
    [ 3  0 -1]
    [ 0 -2  0]
    [ 0  0 -2]
    sage: M.T(3).matrix()
    [ 4  0 -1]
    [ 0 -1  0]
    [ 0  0 -1]
    sage: M.T(5).matrix()
    [ 6  0 -1]
    [ 0  1  0]
    [ 0  0  1]
