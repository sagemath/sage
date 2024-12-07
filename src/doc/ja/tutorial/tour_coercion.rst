.. -*- coding: utf-8 -*-

.. _section-coercion:

================================
ペアレント，型変換および型強制
================================

この節の内容はこれまでと比べるとテクニカルな感じがするかもしれない．
しかし，ペアレントと型強制の意味について理解しておかないと，Sageにおける環その他の代数構造を有効かつ効率的に利用することができないのである．

以下で試みるのは概念の解説であって，それをどうやって実現するかまでは示すことはできない．
実装法に関するチュートリアルは `Sage thematic tutorial <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_ にある．


元
--------

ある環をPythonを使って実装する場合，その第一歩は目的の環の元 ``X`` を表すクラスを作り， ``__add__`` ， ``__sub__`` ， ``__mul__`` のようなダブルアンダースコア メソッド(フックメソッド)によって環の公理を保証する演算を実現することだ．


Pythonは(ダイナミックではあっても)強い型付けがなされる言語なので，最初のうちは環それぞれを一つのPythonクラスで実装すべきだろうと思うかもしれない．
なんと言っても，Pythonは整数については ``<int>`` ，実数については ``<float>`` といった具合に型を一つづつ備えているわけだし．
しかし，このやり方はすぐに行き詰まらないわけにはいかない．
環の種類が無限だからと言って，対応するクラスを無限に実装することはできないからだ．

代りに，群，環，斜体(加除環)，可換環，体，代数系などという順で，普遍的な代数構造の実装を目的とするクラスの階層を作りあげようとする人もいるかもしれない．

しかし，これは明確に異なる環に属する元が同じ型を持ちうることを意味する．
::

    sage: P.<x,y> = GF(3)[]
    sage: Q.<a,b> = GF(4,'z')[]
    sage: type(x)==type(a)
    True


一方，数学的には同等の構造物に対して，(例えば密行列と疎行列に対するように)別々のPythonクラスが異なった形で実装されることもありうる．
::

    sage: P.<a> = PolynomialRing(ZZ)
    sage: Q.<b> = PolynomialRing(ZZ, sparse=True)
    sage: R.<c> = PolynomialRing(ZZ, implementation='NTL')
    sage: type(a); type(b); type(c)
    <class 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
    <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_integral_domain_with_category.element_class'>
    <class 'sage.rings.polynomial.polynomial_integer_dense_ntl.Polynomial_integer_dense_ntl'>


以上から，解決すべき問題は二系統あることが分る．
ある二つの元が同じPythonクラス由来のインスタンスであるとすると，付随する ``__add__`` メソッドによる加算が可能になっているはずだ．
しかし，これら二つが数学的には非常に異なる環に属しているのならば，加算は不能にしておきたい．
一方，数学的に同一の環に属している元に対しては，異なる実装に由来していてもそれらの元の加算は可能であるべきだろう．
異なるPythonクラスに由来する限り，これは簡単に実現できることではない．


これらの問題に対する解は「型強制」(coercion)と呼ばれており，以降で解説する．


しかし，まず肝心なのは，全ての元が自分の帰属先を知っていることだ．
これを可能にするのが ``parent()`` メソッドである:

.. link

::

    sage: a.parent(); b.parent(); c.parent()
    Univariate Polynomial Ring in a over Integer Ring
    Sparse Univariate Polynomial Ring in b over Integer Ring
    Univariate Polynomial Ring in c over Integer Ring (using NTL)


ペアレントとカテゴリー
-------------------------

Pythonが代数構造の元に対応するクラス階層を備えているように，Sageもそれらの元を含む代数構造に対応するクラス群を提供している．
Sageでは元が属する構造物のことを「ペアレント構造」(parent structure)と呼び，基底となるクラスを持つ．
そうしたクラス群は，おおむね数学的概念に沿った形で，集合，環，体などといった順の階層を形成している．

::

    sage: isinstance(QQ,Field)
    True
    sage: isinstance(QQ, Ring)
    True
    sage: isinstance(ZZ,Field)
    False
    sage: isinstance(ZZ, Ring)
    True

代数学では，同じ種類の代数構造を共有する物を，いわゆる「圏」(category)と呼ばれるものに集約して扱う．
Sageのクラス階層と圏の階層構造にはそれなりに類似が見られないでもない．
しかし，Pythonクラスについては圏との類似はあまり強調すべきものでもなさそうだ．
いずれにせよ，数学的な意味における圏はSageでも実装されている:


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

Sageにおけるクラス階層は具体的実装に焦点が当てられている一方，Sageの圏フレームワークではより数学的な構造が重視されている．
圏に対する個々の実装からは独立した，包括的なメソッドとテストを構成することが可能である．


Sageにおけるペアレント構造は，Pythonオブジェクトとして唯一のものであると仮定されている．
例えば，いったんある基底環上の多項式環が生成元と共に作成されると，その結果は実記憶上に配置される:

::

    sage: RR['x','y'] is RR['x','y']
    True


型とペアレント
--------------------

``RingElement`` 型は数学的概念としての「環の元」に完璧に対応しているわけではない．
例えば，正方行列は一つの環に属していると見なしうるにもかかわらず， ``RingElement`` のインスタンスにはならない:


::

    sage: M = Matrix(ZZ,2,2); M
    [0 0]
    [0 0]
    sage: isinstance(M, RingElement)
    False


*ペアレント* が唯一のものであるとしても，同じSageのペアレントに由来する対等な *元* までが同一になるとは限らない．
この辺りはPythonの(全てではないにしても)整数の振舞いとは違っている．

::

    sage: int(1) is int(1) # Pythonのint型
    True
    sage: int(-15) is int(-15)
    False
    sage: 1 is 1           # Sageの整数
    False


重要なのは，異なる環に由来する元は，一般にその型ではなくペアレントによって判別されることである:

::

    sage: a = GF(2)(1)
    sage: b = GF(5)(1)
    sage: type(a) is type(b)
    True
    sage: parent(a)
    Finite Field of size 2
    sage: parent(b)
    Finite Field of size 5

とういうわけで，代数学的な立場からすると **元のペアレントはその型より重要である** ことになる．


型変換と型強制
--------------------------

場合によっては，あるペアレント構造に由来する元を，異なるペアレント構造の元へ変換することができる．
そうした変換は明示的に，あるいは暗黙的に行なうことが可能で，後者を *型強制* (coercion)と呼ぶ．


読者は，例えばC言語における *型変換* (type conversion)と *型強制* (type coercion)の概念をご存知かもしれない．
Sageにも *型変換* と *型強制* の考えは取り込まれている．
しかし，Sageでは主たる対象が型ではなくペアレントになっているので，Cの型変換とSageにおける変換を混同しないよう注意していただきたい．

以下の説明はかなり簡略化されているので，詳しい解説と実装情報についてはSageレファレンスマニュアルの型強制に関する節と `thematic tutorial <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_ を参照されたい．

*異なる* 環に属する元同士の演算実行については，両極をなす二つの立場がある:


* 異なる環はそれぞれが異なる世界を形作っており，何であれ異なる環由来の元同士で和や積を作ることは意味をなさない．
  ``1`` は整数であるのに ``1/2`` が有理数なのだから， ``1 + 1/2`` ですら意味をもちえない．


という立場もあるし

* 環 ``R1`` の元 ``r1`` が何とか他の環 ``R2`` の元と見なしうるなら， ``r1`` と ``R2`` の任意の元に対する全ての算術演算が許される．単位元は全ての体と多くの環に存在し，全て等価と見なしうる．

と考える立場もありうる．


Sageが宗とするのは歩み寄りだ．
``P1`` と ``P2`` がペアレント構造で ``p1`` が ``P1`` の元であるとき， ``p1`` が ``P2`` に帰属するとする解釈をユーザが明示的に求めることがあるかもしれない．
この解釈があらゆる状況で有意であるとは限らないし， ``P1`` の全ての元に対して適用可能とも言えない．
その解釈が意味を持つかどうかはユーザの判断にかかっているのである．
我々はこうした解釈の要求を， **変換** (conversion) と呼ぶことにする:


::

    sage: a = GF(2)(1)
    sage: b = GF(5)(1)
    sage: GF(5)(a) == b
    True
    sage: GF(2)(b) == a
    True


しかし， *暗黙的* (自動的) 変換については，変換が *全面的* かつ *無矛盾* に行ないうる場合にのみ実行される．
こちらで重視されているのは数学的な厳密さである．


そうした暗黙的変換は **型強制** (coercion)と呼ばれる．
型強制が定義できるのならば，結果は型変換と一致しなければならない．
型強制の定義に際して満足されるべき条件は二つある:


#. ``P1`` から ``P2`` への型強制は構造保存写像(すなわち環準同形写像)になっていなければならない．
   ``P1`` の要素が ``P2`` に写像されるだけでは不十分で，その写像は ``P1`` の代数構造を反映している必要がある．

#. 型強制は無矛盾に構成されなければならない．
   ``P3`` を３つ目のペアレント構造として， ``P1`` から ``P2`` への型強制と
   ``P2`` から ``P3`` への型強制を合成すると， ``P1`` から ``P3`` への型強制に一致しなければならない．
   特に ``P1`` から ``P2`` へと ``P2`` から ``P1`` への型強制が存在する場合，この2つの変換を合成すると ``P1`` への恒等写像にならねばならない．


したがって， ``GF(2)`` の全ての元は ``GF(5)`` 上へ変換可能であるにも関わらず，型強制は成立しない．
``GF(2)`` と ``GF(5)`` の間には環準同形写像が存在しないからである．


二つ目の条件 --- 無矛盾性 --- については，いくぶん説明が難しいところがある．
多変数多項式環を例にとって説明してみたい．
実用上，変数名を維持しない型強制はまず使いものにならないはずだ．であれば:


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


変数名を維持する環準同形写像が定義できなければ，型強制も成立しない．
しかし，対象とする環の生成元を生成元リスト上の順序に応じて写像してやれば，型変換の方はまだ定義の可能性が残る:

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

ところが，そうした順序依存の変換は型強制としては満足すべきものにならない．
``ZZ['x','y']`` から ``ZZ['y','x']`` への変数名維持写像と ``ZZ['y','x']`` から ``ZZ['a','b']`` への順序依存写像を合成すると，結果は変数名も順序も保存しない写像となって無矛盾性が破れてしまうからである．


型強制が成立するなら，異なる環に由来する元同士の比較や算術演算の際に利用されるはずである．
これはたしかに便利なのだが，ペアレントの違いを越えた ``==`` 型関係の適用には無理が生じがちなことには注意を要する．
``==`` は *同一の* 環上の元同士の等価関係を表わすが，これは *異なる* 環の元が関わると必ずしも有効なわけではない．
例えば， ``ZZ`` 上の ``1`` と，何か有限体上にあるとした ``1`` は等価であると見なすことができる．
というのは，整数から任意の有限体へは型強制が成り立つからだ．
しかし，一般には二つの異なる有限体環の間に型強制は成立しない．
以下を見ていただきたい:


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


同様にして


.. link

::

    sage: R3(R1.1) == R3.1
    True
    sage: R1.1 == R3.1
    False
    sage: R1.1 != R3.1
    True


さらに無矛盾性の条件から帰結するのは，厳密な環(例えば有理数 ``QQ``)から厳密ではない環(例えば有限精度の実数 ``RR``)への型強制は成立するが，逆方向は成立しないことである．
``QQ`` から ``RR`` への型強制と ``RR`` から ``QQ`` への変換を合成すると ``QQ`` 上の恒等写像になるはずだが，これは不可能である．
と言うのは，有理数の中には，以下で示すように ``RR`` 上で問題なく扱えるものがあるからだ:

::

    sage: RR(1/10^200+1/10^100) == RR(1/10^100)
    True
    sage: 1/10^200+1/10^100 == 1/10^100
    False

型強制が成立しない環 ``P1`` と ``P2`` の二つのペアレント由来の元を比較するとき，基準となるペアレント ``P3`` が選択できて ``P1`` と ``P2`` を ``P3`` へ型強制できる場合がある．
そうした状況では型強制がうまく成立するはずだ．
典型的な例は有理数と整数係数の多項式の和の計算で，結果は有理係数の多項式になる．


::

    sage: P1.<x> = ZZ[]
    sage: p = 2*x+3
    sage: q = 1/2
    sage: parent(p)
    Univariate Polynomial Ring in x over Integer Ring
    sage: parent(p+q)
    Univariate Polynomial Ring in x over Rational Field


この結果は，原則的には ``ZZ['x']`` の有理数体上でも成立する．
しかし，Sageは最も自然に見える *正準* な共通のペアレントを選択しようとする(ここでは ``QQ['x']``)．
共通のペアレント候補が複数あってどれも同じく有望そうな場合，Sageは中の一つをランダムに選択するということは *しない* ．
これは再現性の高い結果を求めるためで，選択の手段については `thematic tutorial
<http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_
に解説がある．


以下に示すのは，共通のペアレントへの型強制が成立しない例である:

::

    sage: R.<x> = QQ[]
    sage: S.<y> = QQ[]
    sage: x+y
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +: 'Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

だめな理由は，Sageが有望そうな候補 ``QQ['x']['y']`` ， ``QQ['y']['x']`` ， ``QQ['x','y']`` あるいは ``QQ['y','x']`` のどれも選択できないことである．
と言うのも，これら4つの相異なる構造はどれも共通なペアレントとして相応しく，基準となるべき選択肢にならないからだ．

