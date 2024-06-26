************
はじめに
************

このチュートリアルは，３〜４時間あればじゅうぶん読み通すことができるはずだ．
HTML版とPDF版のどちらを読んでもいいし、Sageノートブックを経由することもできる(チュートリアル内容をSageから対話的に実行するには，ノートブックで ``Help``,  続けて ``Tutorial`` をクリックする)．

Sageのかなりの部分がPythonを使って実装されているものの，このチュートリアルを読むについてはPythonの予備知識はいらない．
いずれはPythonを勉強したくなるはずだが(とても面白い言語だ)，そんな場合のためにはPythonビギナーズガイド [PyB]_ にリストがある優れた教材がフリーでたくさん用意されている．
とにかく手っ取り早くSageを試してみたいだけなら、このチュートリアルがよい出発点になる．
例えばこんな具合だ:

::

    sage: 2 + 2
    4
    sage: factor(-2007)
    -1 * 3^2 * 223

    sage: A = matrix(4,4, range(16)); A
    [ 0  1  2  3]
    [ 4  5  6  7]
    [ 8  9 10 11]
    [12 13 14 15]

    sage: factor(A.charpoly())
    x^2 * (x^2 - 30*x - 80)

    sage: m = matrix(ZZ,2, range(4))
    sage: m[0,0] = m[0,0] - 3
    sage: m
    [-3  1]
    [ 2  3]

    sage: E = EllipticCurve([1,2,3,4,5]);
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
    over Rational Field
    sage: E.anlist(10)
    [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    sage: E.rank()
    1

    sage: k = 1/(sqrt(3)*I + 3/4 + sqrt(73)*5/9); k
    36/(20*sqrt(73) + 36*I*sqrt(3) + 27)
    sage: N(k)
    0.165495678130644 - 0.0521492082074256*I
    sage: N(k,30)      # 精度は30ビット
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{36}{20 \, \sqrt{73} + 36 i \, \sqrt{3} + 27}

.. _installation:

インストール
==============

まだSageをコンピュータにインストールしていないけれども何かコマンドを実行してはみたいというなら， http://sagecell.sagemath.org 上でオンライン実行してみる手がある．

Sageを自分のコンピュータへインストールする手順については，本家Sageウェブページ [SA]_ のドキュメンテーション部にある "Sage Installation Guide"を見てほしい．
ここではいくつかコメントしておくだけにしよう．

#. Sageのダウンロード用ファイルは「バッテリー込み」である．
   つまり、SageはPython, IPython, PARI, GAP,  Singular, Maxima, NTL, GMPなどを援用して動作するが，これらは全てSageの配布ファイルに含まれているので別途インストールする必要はないということだ．
   ただし、MacaulayやKASHなど一部の機能を利用するには関連するオプショナルなSageパッケージをインストールするか，少なくとも使うコンピュータに関連プログラム群がインストール済みでなくてはならない(利用できるオプショナル・パッケージの一覧を見るには, ``sage -optional`` を実行するか，あるいはSageウェブサイトの "Download"ページを見るとよい)．

#. コンパイル済みのバイナリ版Sage(Sageサイトにある)は，ソースコードより簡単かつ速やかにインストールすることができる．
   ファイルを入手したら展開して ``sage`` コマンドを実行するだけで出来上がりだ．

#. SageTeXパッケージを使いたいのならば(SageTeXはSageの処理結果をLaTeX文書に埋め込み可能にしてくれる)，使用すべきTeXディストリビューションをSageTeXに教えてやる必要がある．
   設定法については， `Sage installation guide <http://doc.sagemath.org/html/en/>`_ 中の "Make SageTeX known to TeX" を参照してほしい(ローカルシステム上の `ここ <../../en/installation/index.html>`_ にもインストールガイドがある)．
   手順はごく簡単で，環境変数を一つ設定するか，あるいはTeX配下のディレクトリにファイルを1個コピーしてやるだけである．


SageTeXの利用に関する解説は
``$SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/`` にある．
``$SAGE_ROOT`` はSageがインストールされているディレクトリで，例えば ``/opt/sage-9.6`` などとなっているはずだ．


Sageの使いかた
================

Sageを使うには以下のようなやり方がある．

- **ノートブック グラフィカル インターフェイス:**  ``sage -n jupyter`` を実行する.
  `Jupyter documentation on-line <https://jupyter-notebook.readthedocs.io/en/latest/notebook.html>`_ を読む.

- **対話的コマンドライン:** :ref:`chapter-interactive_shell` 節を参照．

- **プログラム作成:** Sage上でインタープリタおよびコンパイラを経由してプログラムを書く(:ref:`section-loadattach` 節と :ref:`section-compile` 節を参照)．さらに

- **スクリプト作成:** Sageライブラリを利用するスタンドアロンPythonスクリプトを書く(:ref:`section-standalone` 節を参照).


Sageの長期目標
=======================

- **有用性**: Sageが想定しているユーザは，数学を学ぶ学生(高校生から大学学部生まで)と教師、そして数学の専門家である．
  代数、幾何、数論、解析学、数値解析などの数学諸分野には，種々の概念や量が現われてくる．
  Sageの狙いは、ユーザが数学上の概念や諸量の性質を探ったり，それらの働きを体験する手助けになるようなソフトウェアを提供することである．
  Sageを使えば，各種の数学的な実験を容易に対話的に実行することができる．

- **高速性:** 動作が高速である．
  Sageは GMP, PARI, GAP, NTLなど高度に最適化された完成度の高いソフトウェアを援用しており，多くの場合きわめて高速に演算が実行される．

- **フリーかつオープンソース:** ソースコードは自由に入手可能で，可読性が高くなければならない．
  そうすればユーザはSageが行なう処理の詳細を理解することができるし，拡張も容易になる．
  数学者であれば，定理を深く理解するために証明をていねいに読むか，少なくとも証明の流れ程度は追っておくはずである．
  計算システムのユーザも同じことで，演算処理がどのように実行されるのかソースコードを読んで把握できるようであってほしい．
  論文発表する仕事の計算にSageを使っておけば，論文の読者も確実にSageとその全ソースコードを自由に利用できることになる．
  Sageでは，仕事に使ったバージョンを保存しておいて再配布することすら許されているのだ．

- **コンパイルが容易:** Sageは，Linux， OSXあるいはWindowsのユーザがソースコードから容易にコンパイル・ビルドできるようでなくてはならない．
  これによりユーザはSageシステムを柔軟に修正することができる．

- **協調性:** Sageは，PARI， GAP， Singular， Maxima， KASH， Magma， Maple，さらにMathematicaなど多くのコンピュータ代数システムとの頑健なインターフェイスを提供する．
  Sageの狙いは、既存の数学ソフトウェアとの統合と拡張である．

- **豊富な関連文書:** チュートリアル，プログラミングガイド，レファレンスマニュアル，ハウツー類が揃っている．
  これには多数の具体例と数学的背景知識の解説も含まれる．

- **拡張性:** 新しいデータ型をゼロから定義したり，既存のデータ型を利用して作り出すことができる．
  さまざまな言語で書いたプログラムをシステムに組み込んで利用することも可能だ．

- **ユーザーフレンドリー**: ユーザは使用するオブジェクトにどんな属性や機能が組込まれているかを簡単に把握し，さらに関連文書やソースコードなども容易に閲覧できなくてはならない．
  高度のユーザーサポートも提供される．
