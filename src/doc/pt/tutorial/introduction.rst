**********
Introdução
**********

Este tutorial leva no máximo de 3 a 4 horas para ser percorrido. Você
pode lê-lo em versão HTML ou PDF, ou a partir do Notebook Sage (clique
em ``Help``, então clique em ``Tutorial`` para percorrer o tutorial de
forma interativa).

Embora grande parte do Sage seja implementado em Python, nenhum
conhecimento de Python é necessário para a leitura deste tutorial.
Você vai querer aprender Python (uma linguagem muito divertida!) em
algum momento, e existem diversas opções gratuitas disponíveis para
isso, Guia do Iniciante em Python [PyB]_ (em inglês) lista muitas opções.
Se você quiser experimentar o Sage rapidamente, este tutorial é o lugar
certo para começar. Por exemplo:

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
    sage: N(k,30)      # 30 "bits"
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{36}{20 \, \sqrt{73} + 36 i \, \sqrt{3} + 27}

.. _installation:

Instalação
==========

Se você não tem o Sage instalado em um computador e quer apenas
experimentar alguns comandos, use o Sage através do site
http://sagecell.sagemath.org.

Veja o guia de instalação do Sage na seção de documentação na página
principal do Sage [SA]_ para instruções de como instalar o Sage no seu
computador. Aqui faremos apenas alguns comentários.

#. O arquivo para instalação do Sage vem com "baterias incluídas". Em
   outras palavras, embora o Sage use o Python, IPython, PARI, GAP,
   Singular, Maxima, NTL, GMP, e uma série de outros programas, você
   não precisa instalá-los separadamente pois eles estão incluídos no
   Sage. Todavia, para usar alguns recursos, por exemplo, o Macaulay
   ou o KASH, você precisa
   ter os programas necessários já instalados no seu computador.

#. A versão pré-compilada do Sage (disponível na página do Sage na
   internet) pode ser mais fácil e rápida para instalar do que a
   versão obtida compilando o código fonte.

#. Se você quiser usar o pacote SageTeX (que permite inserir
   cálculos do Sage em um arquivo LaTeX), você deve tornar
   o SageTex disponível para a sua distribuição TeX. Para fazer isso,
   consulte a seção "Make SageTex known to TeX" no `Sage installation
   guide <http://doc.sagemath.org/html/en/>`_. O procedimento é bem
   simples; você precisa apenas definir algumas variáveis no seu
   sistema ou copiar um arquivo para um diretório onde o TeX poderá
   encontrá-lo.

   A documentação para usar o SageTex está disponível em
   ``$SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/``, onde
   ``$SAGE_ROOT`` refere-se ao diretório onde você instalou o Sage
   -- por exemplo, ``/opt/sage-9.6``.

Formas de usar o Sage
=====================

Você pode usar o Sage de diversas formas.


-  **Interface gráfica Notebook:** inicie ``sage -n jupyter``; leia
   `Jupyter documentation on-line <https://jupyter-notebook.readthedocs.io/en/latest/notebook.html>`_,

-  **Linha de comando interativa:** veja
   :ref:`chapter-interactive_shell`,

-  **Programas:** escrevendo programas interpretados e compilados em
   Sage (veja :ref:`section-loadattach` e :ref:`section-compile`), e

-  **Scripts:** escrevendo scripts em Python que usam a biblioteca do
   Sage (veja :ref:`section-standalone`).


Objetivos do Sage a longo prazo
===============================

-  **Útil**: O público alvo do Sage são estudantes de matemática
   (desde o ensino médio até a pós-graduação), professores, e
   pesquisadores em matemática. O objetivo é fornecer um software que
   possa ser usado para explorar e experimentar construções matemáticas
   em álgebra, geometria, teoria de números, cálculo, computação
   numérica, etc. O Sage torna mais fácil a experimentação com objetos
   matemáticos de forma interativa.

-  **Eficiente:** Ser rápido. O Sage usa software bastante otimizado
   como o GMP, PARI, GAP, e NTL, e portanto é muito rápido em certas
   operações.

-  **Gratuito e de código aberto:** O código fonte deve ser amplamente
   disponível e legível, de modo que os usuários possam entender o que
   o software realmente faz e possam facilmente estendê-lo. Da mesma
   forma que matemáticos ganham entendimento sobre um teorema lendo
   cuidadosamente a sua demonstração, as pessoas que fazem cálculos
   deveriam poder entender como os cálculos são feitos lendo o código
   fonte e seus comentários. Se você usar o Sage para fazer cálculos em
   um artigo que seja publicado, você pode ter certeza que os leitores
   sempre terão livre acesso ao Sage e seu código fonte, e você tem até
   mesmo permissão para arquivar e redistribuir a versão do Sage que
   você utilizou.

-  **Fácil de compilar:** O Sage deve ser fácil de compilar a partir
   do código fonte para usuários de Linux, OS X e Windows. Isso
   fornece mais flexibilidade para os usuários modificarem o sistema.

-  **Cooperação:** Fornecer uma interface robusta para outros sistemas
   computacionais, incluindo PARI, GAP, Singular, Maxima, KASH, Magma,
   Maple e Mathematica. O Sage foi concebido para unificar e estender
   outros softwares de matemática existentes.

-  **Bem documentado:** Tutorial, guia de programação, manual de
   referência, e how-to, com inúmeros exemplos e discussão sobre
   conceitos matemáticos relacionados.

-  **Estensível:** Ser capaz de definir novos tipos de dados ou
   derivá-los a partir dos tipos de dados existentes, e usar programas
   escritos em diversas outras linguagens.

-  **Fácil de usar:** Deve ser fácil entender quais recursos estão
   disponíveis para um determinado objeto e consultar a documentação e
   o código fonte.
