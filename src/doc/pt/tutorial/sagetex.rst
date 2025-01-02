.. _sec-sagetex:

****************
Usando o SageTeX
****************

O pacote SageTeX permite que você insira resultados de cálculos feitos
com o Sage em um documento LaTeX. Esse pacote já vem com o Sage. Para
usá-lo, você precisa "instalá-lo" em seu sistema LaTeX local; aqui
instalar significa copiar um simples arquivo. Veja :ref:`installation`
neste tutorial e a seção "Make SageTeX known to TeX" do `Guia de
instalação do Sage <http://doc.sagemath.org/html/en/installation/index.html>`_
(em inglês).

Aqui vai um breve exemplo de como usar o SageTeX. A documentação
completa pode ser encontrada em
:file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex`, onde
``SAGE_ROOT`` é o diretório onde se encontra a sua instalação. Esse
diretório contém a documentação, um arquivo de exemplo, e alguns
scripts em Python possivelmente úteis.

Para ver como o SageTeX funciona, siga as instruções para instalar o
SageTeX (em :ref:`installation`) e copie o seguinte texto em um
arquivo chamado ``st_example.tex``, por exemplo.

.. warning::

  O texto abaixo vai apresentar diversos erros sobre "unknown control
  sequences" se você está visualizando isto na ajuda "live". Use a
  versão estática para ver o texto corretamente.

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

Execute o LaTeX em ``st_example.tex`` da forma usual. Note que o LaTeX
vai reclamar sobre algumas coisas, entre elas::

    Package sagetex Warning: Graphics file
    sage-plots-for-st_example.tex/plot-0.eps on page 1 does not exist. Plot
    command is on input line 25.

    Package sagetex Warning: There were undefined Sage formulas and/or
    plots. Run Sage on st_example.sage, and then run LaTeX on
    st_example.tex again.

Observe que, além dos arquivos usuais produzidos pelo LaTeX, existe um
arquivo chamado ``st_example.sage``. Esse é um script em Sage
produzido quando você executa o LaTeX em ``st_example.tex``. A
mensagem de alerta pede para você executar o LaTeX em
``st_example.sage``, então siga essa sugestão e faça isso. Você vai
receber uma mensagem para executar o LaTeX em ``st_example.tex``
novamente, mas antes que você faça isso, observe que um novo arquivo
foi criado: ``st_example.sout``. Esse arquivo contém os resultados dos
cálculos feitos pelo Sage, em um formato que o LaTeX pode usar para
inserir em seu texto. Um novo diretório contendo um arquivo EPS do seu
gráfico também foi criado. Execute o LaTeX novamente e você vai ver
que tudo que foi calculado, incluindo os gráficos, foi incluído em seu
documento.

As funções (macros em inglês) utilizadas acima devem ser fáceis de
entender. Um ambiente ``sageblock`` insere código "verbatim"
(exatamente como é digitado) e também executa o código quando você
executa o Sage. Quando você insere ``\sage{foo}``, é incluído em seu
documento o resultado que você obteria executando ``latex(foo)`` no
Sage. Comandos para fazer gráficos são um pouco mais complicados, mas
em sua forma mais simples, ``\sageplot{foo}`` insere a imagem que você
obtêm usando ``foo.save('filename.eps')``.

Em geral, a rotina é a seguinte:

    - execute o LaTeX no seu arquivo .tex;
    - execute o Sage no arquivo .sage que foi gerado;
    - execute o LaTeX novamente.

Você pode omitir a execução do Sage desde que você não tenha alterado
os comandos em Sage em seu documento.

Há muito mais sobre o SageTeX, e como tanto o Sage como o LaTeX são
ferramentas complexas e poderosas, é uma boa idéia ler a documentação
para o SageTeX que se encontra em
:file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex`.
