.. _sec-sagetex:

===============
 SageTeXを使う
===============

SageTeXパッケージを使うと，Sageによる処理結果をLaTeX文書に埋め込むことができるようになる．
利用するためには，まずSageTeXを「インストール」しておかなければならない(:ref:`sec-sagetex_install` 節を参照)．


具体例
=======

ここでは，ごく簡単な例題を通してSageTeXの利用手順を紹介する．
完全な解説ドキュメントと例題ファイルは，ディレクトリ :file:`SAGE_ROOT/venv/share/doc/sagetex` に置いてある．
:file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex` にあるPythonスクリプトは何か役に立つ場面があるはずだ．
以上の ``SAGE_ROOT`` は，Sageをインストールしたディレクトリである．


SageTeXの動作を体験するために，まずSageTeXのインストール手続き(:ref:`sec-sagetex_install` 節)を実行し，
以下のテキストを ``st_example.tex`` などという名前で保存しておいてほしい:


.. warning::

  このテキストを "live"ヘルプから表示すると命令未定義エラーになる．
  正常に表示するためには "static"ヘルプで表示すること．


.. code-block:: latex

    \documentclass{article}
    \usepackage{sagetex}

    \begin{document}

    Using Sage\TeX, one can use Sage to compute things and put them into
    your \LaTeX{} document. For example, there are
    $\sage{number_of_partitions(1269)}$ integer partitions of $1269$.
    You don't need to compute the number yourself, or even cut and paste
    it from somewhere.

    Here's some Sage code:

    \begin{sageblock}
        f(x) = exp(x) * sin(2*x)
    \end{sageblock}

    The second derivative of $f$ is

    \[
      \frac{\mathrm{d}^{2}}{\mathrm{d}x^{2}} \sage{f(x)} =
      \sage{diff(f, x, 2)(x)}.
    \]

    Here's a plot of $f$ from $-1$ to $1$:

    \sageplot{plot(f, -1, 1)}

    \end{document}

この ``st_example.tex`` をいつも通りにLaTeXで処理する．
するとLaTeXは以下のような文句をつけてくるだろう:


::

    Package sagetex Warning: Graphics file
    sage-plots-for-st_example.tex/plot-0.eps on page 1 does not exist. Plot
    command is on input line 25.

    Package sagetex Warning: There were undefined Sage formulas and/or
    plots. Run Sage on st_example.sagetex.sage, and then run LaTeX on
    st_example.tex again.

注目してほしいのは，LaTeXが通常の処理で生成するファイル群に加えて， ``st_example.sage`` というファイルが出来ていることだ．
これは ``st_example.tex`` の処理畤に生成されたSageスクリプトで，上で見たLaTeX処理時のメッセージは，この ``st_example.sage`` をSageで実行せよという内容である．
その通り実行すると ``st_example.tex`` を再びLaTeXで処理せよと告げられるが，その前に新しいファイル ``st_example.sout`` が生成されていることに注意．
このファイルにはSageの演算結果がLaTeXテキストに挿入して利用可能な形式で保存されている．
プロット画像のEPSファイルを含むディレクトリも新規作成されている．
ここでLaTeX処理を実行すると，Sageの演算結果とプロットの全てがLaTeX文書に収められることになる．


上の処理に用いられた各マクロの内容はごく簡単に理解できる．
``sageblock`` 環境はSageコードを入力通りに組版し，ユーザーがSageを動かすとそのコードを実行する．
``\sage{foo}`` とすると， Sage上で ``latex(foo)`` を実行したのと同じ結果がLaTeX文書に挿入される．
プロット命令はやや複雑だが，もっとも単純な場合である ``\sageplot{foo}`` は ``foo.save('filename.eps')`` を実行して得られた画像を文書へ挿入する役割を果たす．


要するに，必要な作業は以下の三段階になる:

- LaTeXで .texファイルを処理
- 生成された .sageファイルをSageで実行
- LaTeXで .texファイルを再処理


作業中にLaTeX文書内のSageコマンドを変更しない場合，Sageによる処理は省略することができる．


SageTeXは到底以上で語り尽せるものでなく，SageとLaTeXは共に複雑で強力なツールだ．
:file:`SAGE_ROOT/venv/share/doc/sagetex` にあるSageTeXのドキュメントを読むことを強くお勧めする．


.. _sec-sagetex_install:


TeXにSageTeXの存在を教える
===========================

Sageはおおむね自己完結的なシステムなのだが，正しく機能するために外部ツールの介入を要する部分があることも確かだ．
SageTeXもそうした部分の一つである．


SageTeXパッケージを使えばSageによる演算やプロットをLaTeX文書に埋め込むことが可能になる．
SageTeXはデフォルトでSageにインストールされるが，LaTeX文書で利用する前に，運用しているTeXシステムへSageTeXの存在を教えておかねばならない．


鍵になるのは， TeXが ``sagetex.sty`` を発見できるかどうかである．
この ``sagetex.sty`` は， ``SAGE_ROOT`` をSageがビルトあるいはインストールされたディレクトリとすると，
:file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/` に置かれているはずだ．
TeXが ``sagetex.sty`` を読めるようにしてやらなければ，SageTeXも動作できないのである．
これを実現するには何通りかのやり方がある．


- 第一の，かつ一番簡単な方法は， ``sagetex.sty`` を作成すべきLaTeX文書と同じディレクトリ内にコピーしておくことである．
  TeXは組版処理の際に現ディレクトリを必ずサーチするから，この方法は常に有効だ．

  ただし，このやり方には二つのちょっとした問題点がある．
  一つ目は，このやり方では使用しているシステムが重複した ``sagetex.sty`` だらけになってしまうこと．
  二つ目の，もっと厄介な問題は，この状態でSageが更新されてSageTeXも新しいバージョンになった場合，SageTeXを構成するPythonコードやLaTeXコードとの食い違いが生じて実行時にエラーが発生しかねない点である．


- 第二の方法は，環境変数 ``TEXINPUTS`` を利用することである．
  bashシェルを使っているなら

  ::

      export TEXINPUTS="SAGE_ROOT/venv/share/texmf//:"

  と実行すればよい．ただし ``SAGE_ROOT`` はSageのインストール先ディレクトリである．
  上の実行例では，行末にスラッシュ2個とコロンを付け忘れないでいただきたい．
  実行後は，TeXと関連ツールがSageTeXスタイルファイルを見つけられるようになる．
  上のコマンド行を ``.bashrc`` に付加して保存しておけば設定を永続させることができる．
  bash以外のシェルを使っている場合， ``TEXINPUTS`` 変数を設定するためのコマンドも異なる可能性がある．
  設定法については，自分の使っているシェルのドキュメントを参照のこと．

  この方法にも瑕はある．
  ユーザがTeXShopやKile，あるいはEmacs/AucTeXなどを使っている場合，必ずしも環境変数を認識してくれるとは限らないのである．
  これらのアプリケーションが常にシェル環境を通してLaTeXを起動するわけではないからだ．

  インストール済みのSageを移動したり，新バージョンを旧版とは違う場所にインストールした場合，
  先に紹介したコマンドも新しい ``SAGE_ROOT`` を反映させるように変更する必要がある．


- TeXに ``sagetex.sty`` の在処を教える第三の(かつ最善の)方法は，このスタイルファイルを自分のホームディレクトリのどこか都合のよい所にコピーしておくことだ．
  TeXディストリビューションの多くは，パッケージを求めてホームディレクトリにある ``texmf`` ディレクトリを自動的に探索するようになっている．
  このディレクトリを正確に特定するには，コマンド

  ::

      kpsewhich -var-value=TEXMFHOME

  を実行する．すると ``/home/drake/texmf`` や ``/Users/drake/Library/texmf`` などと表示されるはずだから， :file:`SAGE_ROOT/venv/share/texmf/` 内の ``tex/`` ディレクトリをホームディレクトリの ``texmf`` にコピーするには

  ::

      cp -R SAGE_ROOT/venv/share/texmf/tex TEXMFHOME

  などとする．
  もちろん， ``SAGE_ROOT`` を実際にSageをインストールしたディレクトリとするのはこれまでと同じことで， ``TEXMFHOME`` は上で見た ``kpsewhich`` コマンドの結果で置き換える．

  SageをアップグレードしたらSageTeXがうまく動かなくなったという場合は，上記の手順をもう一度繰り返すだけでSageTeXのSageとTeX関連部分が同期する．


.. _sagetex_installation_multiuser:

- 複数ユーザに対応するシステムでは，以上の手続きを変更して ``sagetex.sty`` を公開運用中のTeXディレクトリにコピーすればよい．
  おそらく一番賢いコピー先は ``TEXMFHOME`` ディレクトリではなく，コマンド

  ::

      kpsewhich -var-value=TEXMFLOCAL


  の実行結果に従うことだろう．出力は ``/usr/local/share/texmf`` のようになるはずで， 上と同じように ``tex`` ディレクトリを ``TEXMFLOCAL`` ディレクトリ内にコピーする．
  ついでTeXのパッケージデータベースを更新しなければならないが，これは簡単で，ルート権限で

  ::

      texhash TEXMFLOCAL


  と実行すればよい．ただし ``TEXMFLOCAL`` を現実に合わせて変更するのは先と同じだ．
  これでシステムの全ユーザはSageTeXパッケージへアクセス可能になり，Sageが利用できればSageTeXも使えるようになる．

.. warning::

  肝心なのは，LaTeXが組版処理時に使う ``sagetex.sty`` ファイルと，Sageが援用するSageTeXのバージョンが一致していることである．
  Sageを更新したら，あちこちに散らばった古いバージョンの ``sagetex.sty`` を面倒でも全て削除してやらなければいけない．

  SageTeX関連ファイルをホームディレクトリの ``texmf`` ディレクトリ内にコピーしてしまうこと(先に紹介した第三の方法)をお勧めするのは，この面倒があるからである．
  第三の方法にしておけば，Sage更新後もSageTeXを正常に動作させるために必要な作業はディレクトリを一つコピーするだけになる．


SageTeXドキュメント
---------------------

厳密にはSageのインストール一式には含まれないものの，ここで
SageTeXのドキュメントが :file:`SAGE_ROOT/venv/share/doc/sagetex/sagetex.pdf` に配置されていることに触れておきたい．
同じディレクトリには例題ファイルと，これをLaTeXとSageTeXによってすでに組版処理した結果も用意されている(``example.tex`` と ``example.pdf`` を参照)．
これらのファイルは `SageTeX ページ <https://github.com/sagemath/sagetex>`_ からダンロードすることもできる．


SageTeXとTeXLive
-------------------

混乱を招きかねない問題点の一つとして，人気あるTeXディストリビューション
`TeXLive 2009 <http://www.tug.org/texlive/>`_ にSageTeXが含まれている現実があげられる．
これは有り難い感じがするかもしれないが，SageTeXに関して重要なのはSageとLaTeXの各要素が同期していることだ．
SageとSageTeXは共に頻繁にアップデートされるがTeXLiveはそうではないから，その「同期」のところで問題が生じる．
この文の執筆時点(2013年3月)では，多くのLinuxディストリビューションが新しいTeXLiveリリースに移行しつつある．
しかし2009リリースもしぶとく生き残っていて，実はこれがSageTeXに関するバグレポートの主要な発生源になっているのだ．

このため *強く推奨* させていただきたいのは，SageTeXのLaTeX関連部分は以上で説明したやり方で常にSageからインストールすることである．
上記の手順に従えば，SageTeXのSageおよびLaTeX対応部分の互換性が保証されるから，動作も正常に保たれる．
SageTeXのLaTeX対応部分をTeXLiveから援用することはサポート対象外になる．


