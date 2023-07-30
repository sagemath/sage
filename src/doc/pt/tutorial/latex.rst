*********************************
Sage, LaTeX e Companheiros
*********************************

AUTOR:  Rob Beezer (2010-05-23)

O Sage e o dialeto LaTeX do TeX tem um relacionamento sinergético
intenso. Esta seção tem como objetivo introduzir as diversas formas de
interação entre eles, começado pelas mais básicas e indo até as menos
usuais. (Logo você pode não querer ler esta seção inteira na sua
primeira passagem por este tutorial.)

Panorama Geral
==============

Pode ser mais fácil entender os vários usos do LaTeX com um panorama
geral sobre os três principais métodos usados pelo Sage.

#. Todo objeto no Sage possui uma representação em LaTeX. Você
   pode acessar essa representação executando, no Notebook ou na
   linha de comando do Sage, ``latex(foo)`` where ``foo`` é algum
   objeto no Sage. O resultado é uma string que deve fornecer uma
   representação razoável de ``foo`` no modo matemático em LaTeX
   (por exemplo, quando cercado por um par de símbolos $). Alguns
   exemplos disso seguem abaixo.

   Dessa forma, o Sage pode ser usado efetivamente para construir
   partes de um documento LaTeX: crie ou calcule um objeto no
   Sage, imprima ``latex()`` do objeto e copie-e-cole o resultado
   no seu documento.

#. A interface Notebook é configurada para usar o `MathJax
   <http://www.mathjax.org/>`_ para representar
   fórmulas matemáticas de forma clara em um web navegador. O MathJax é
   uma coleção de rotinas em JavaScript e fontes associadas.
   Tipicamente esses fontes ficam armazenadas em um servidor e são
   enviadas para o navegador juntamente com a página onde elas estão
   sendo usadas. No caso do Sage, o Notebook está sempre conectado a
   um servidor usado para executar os comando do Sage, e esse servidor
   também fornece as fontes do MathJax necessárias. Logo não é
   necessário configurar nada mais para ter formulas matemáticas
   representadas no seu navegador quando você usa o Notebook do Sage.

   O MathJax é implementado para representar um subconjunto grande,
   mas não completo, do TeX. Ele não suporta objetos como, por
   exemplo, tabelas complicadas e seções, e é focado para
   representar acuradamente pequenas fórmulas em TeX. A
   representação automática de fórmulas matemáticas no Notebook é
   obtida convertendo a representação ``latex()`` de um objeto
   (como descrito acima) em uma forma de HTML mais adequada ao
   MathJax.

   Como o MathJax usa as suas próprias fontes de tamanho variável,
   ele é superior a outros métodos que convertem equações, ou
   outros pequenos trechos de TeX, em imagens estáticas.

#. Na linha de comando do Sage, ou no Notebook quando o código em
   LaTeX é complicado demais para o MathJax processar, uma
   instalação local do LaTeX pode ser usada. O Sage inclui quase
   tudo que você precisa para compilar e usar o Sage, mas uma
   exceção significativa é o TeX. Então nessas situações você
   precisa ter o TeX instalado, juntamente com algumas ferramentas
   de conversão, para usar os recursos completos.

Aqui nós demonstramos alguns usos básicos da função ``latex()``. ::

    sage: var('z')
    z
    sage: latex(z^12)
    z^{12}
    sage: latex(integrate(z^4, z))
    \frac{1}{5} \, z^{5}
    sage: latex('a string')
    \text{\texttt{a{ }string}}
    sage: latex(QQ)
    \Bold{Q}
    sage: latex(matrix(QQ, 2, 3, [[2,4,6],[-1,-1,-1]]))
    \left(\begin{array}{rrr}
    2 & 4 & 6 \\
    -1 & -1 & -1
    \end{array}\right)

A funcionalidade básica do MathJax é em sua maior parte automática no
Notebook, mas nós podemos demonstrar esse suporte parcialmente com a
classe ``MathJax``. A função ``eval`` dessa classe converte um objeto
do Sage em sua representação LaTeX e adiciona HTML que por sua vez
evoca a classe "matemática" do CSS, a qual então emprega o MathJax. ::

    sage: from sage.misc.html import MathJax
    sage: js = MathJax()
    sage: var('z')
    z
    sage: js(z^12)
    <html>\[z^{12}\]</html>
    sage: js(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}\]</html>
    sage: js(ZZ[x])
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}[x]\]</html>
    sage: js(integrate(z^4, z))
    <html>\[\frac{1}{5} \, z^{5}\]</html>

Uso Básico
==========

Como indicado acima, a forma mais simples de explorar o suporte do
Sage para o LaTeX é usando a função ``latex()`` para criar código
LaTeX para representar objetos matemáticos. Essas strings podem então
ser incorporadas em documentos LaTeX. Isso funciona igualmente no
Notebook ou na linha de comando do Sage.

No outro extremo está o comando ``view()``. Na linha de comando do
Sage o comando ``view(foo)`` irá criar a representação em LaTeX de
``foo``, incorporar isso em um documento simples em LaTeX, e então
processar o documento usando o LaTeX em seu sistema. Por fim, o
visualizador apropriado será aberto para apresentar o documento
gerado. Qual versão do TeX é usada, e portanto as opções para a saída
e visualizador, podem ser personalizados (veja
:ref:`sec-custom-processing`).

No Notebook, o comando ``view(foo)`` cria uma combinação apropriada de
HTML e CSS para que o MathJax mostre a representação em LaTeX na folha
de trabalho. Para o usuário, ele simplesmente cria uma versão
cuidadosamente formatada do resultado, distinta da saída padrão em
modo texto do Sage. Nem todo objeto no Sage possui uma representação
em LaTeX adequada às capacidades limitadas do MathJax. Nesses casos, a
interpretação pelo MathJax pode ser deixada de lado, e com isso o LaTeX
do sistema é chamado, e o resultado dessa chamada é convertido em uma
imagem que é inserida na folha de trabalho. Como alterar e controlar
esse processo é discutido abaixo na seção
:ref:`sec-custom-generation`.

O comando interno ``pretty_print()`` ilustra a conversão de objetos do
Sage para HTML que emprega o MathJax no Notebook. ::

    sage: pretty_print(x^12)
    x^12
    sage: pretty_print(integrate(sin(x), x))
    -cos(x)

O Notebook tem outros dois recursos para empregar o TeX. O primeiro é
o botão "Typeset" bem acima da primeira célula da folha de trabalho, à
direita dos quatro menus de opções. Quando selecionado, o resultado de
qualquer cálculo vai ser interpretado pelo MathJax. Note que esse
efeito não é retroativo -- células calculadas anteriormente precisam
ser recalculadas para ter o resultado representado pelo MathJax.
Essencialmente, selecionar o botão "Typeset" é equivalente a aplicar o
comando ``view()`` ao resultado de cada célula.

Um segundo recurso disponível no Notebook é possibilidade de inserir
código TeX para fazer anotações na folha de trabalho. Quando o cursos
esta posicionado entre células de modo que uma barra azul fica
visível, então shift+clique irá abrir um mini processador de texto,
TinyMCE. Isso permite digitar texto, usando um editor WSISYG para
criar HTML e CSS. Logo é possível inserir texto formatado para
complementar a folha de trabalho. Todavia, texto entre símbolos $, ou
$$, é interpretado pelo MathJax como "inline" ou "display math"
espectivamente.

.. _sec-custom-generation:

Personalizando a Criação de Código LaTeX
========================================

Exitem várias formas de personalizar o código LaTeX gerado pelo
comando ``latex()``. No Notebook e na linha de comando existe um
objeto pré-definido chamado ``latex`` que possui diversos métodos, os
quais você pode listar digitando ``latex.``, seguido da tecla tab
(note a presença do ponto).

Um bom exemplo é o método ``latex.matrix_delimiters``. Ele pode ser
usado para alterar a notação de matrizes -- parênteses grandes,
colchetes, barras verticais. Nenhuma noção de estilo é enfatizada,
você pode configurar como desejado. Observe como as barras invertidas
usadas em LaTeX requerem uma barra adicional para que elas não sejam
interpretadas pelo Python como um comando (ou seja, sejam implementadas
simplesmente como parte de uma string. ::

    sage: A = matrix(ZZ, 2, 2, range(4))
    sage: latex(A)
    \left(\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right)
    sage: latex.matrix_delimiters(left='[', right=']')
    sage: latex(A)
    \left[\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right]
    sage: latex.matrix_delimiters(left='\\{', right='\\}')
    sage: latex(A)
    \left\{\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right\}

O método ``latex.vector_delimiters`` funciona de forma similar.

A forma como anéis e corpos comuns podem ser representados pode ser
controlada pelo método ``latex.blackboard_bold``. Esses conjuntos são
representados por padrão em negrito, mas podem opcionalmente ser
escritos em letras duplas como é comum em trabalhos escritos. Isso é
obtido redefinindo a macro ``\Bold{}`` que faz parte do Sage. ::

    sage: latex(QQ)
    \Bold{Q}
    sage: from sage.misc.html import MathJax
    sage: js = MathJax()
    sage: js(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}\]</html>

    sage: latex.blackboard_bold(True)
    sage: js(QQ)
    <html>\[\newcommand{\Bold}[1]{\mathbb{#1}}\Bold{Q}\]</html>
    sage: latex.blackboard_bold(False)

É possível aproveitar os recursos do TeX adicionando novas funções
(macros em inglês) e novos pacotes. Primeiro, funções individuais podem
ser adicionadas para serem usadas quando o MathJax interpreta pequenos
trechos de códigos TeX no Notebook. ::

    sage: latex.extra_macros()
    ''
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: latex.extra_macros()
    '\\newcommand{\\foo}{bar}'
    sage: var('x y')
    (x, y)
    sage: latex(x+y)
    x + y
    sage: from sage.misc.html import MathJax
    sage: js = MathJax()
    sage: js(x+y)
    <html>\[\newcommand{\foo}{bar}x + y\]</html>

Macros adicionais usadas dessa forma serão também usadas eventualmente
se a versão do TeX no seu sistema for usada para lidar com algo muito
complicado para o MathJax. O comando ``latex_extra_preamble`` é usado
para construir o preambulo de um documento completo em LaTeX.
Ilustramos a seguir como fazer isso. Novamente note a necessidade de
barras invertidas duplas nas strings do Python. ::


    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: from sage.misc.latex import latex_extra_preamble
    sage: print(latex_extra_preamble())
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: print(latex_extra_preamble())
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    \newcommand{\foo}{bar}

Novamente, para expressões grandes ou mais complicadas do LaTeX, é
possível adicionar pacotes (ou qualquer outra coisa) ao preambulo do
arquivo LaTeX. Qualquer coisa pode ser incorporada no preambulo com o
comando ``latex.add_to_preamble``, e o comando mais especializado
``latex.add_package_to_preamble_if_available`` irá primeiro verificar
se certo pacote está realmente disponível antes de adicioná-lo ao
preambulo

Agora adicionamos o pacote geometry ao preambulo e usamos ele para
definir o tamanho da região na página que o TeX vai usar
(efetivamente definido as margens). Novamente, observe a necessidade
de barras duplas nas strings do Python. ::


    sage: from sage.misc.latex import latex_extra_preamble
    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: latex.add_to_preamble('\\usepackage{geometry}')
    sage: latex.add_to_preamble('\\geometry{letterpaper,total={8in,10in}}')
    sage: latex.extra_preamble()
    '\\usepackage{geometry}\\geometry{letterpaper,total={8in,10in}}'
    sage: print(latex_extra_preamble())
    \usepackage{geometry}\geometry{letterpaper,total={8in,10in}}
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}

Um pacote pode ser adicionado juntamente com a verificação de sua
existência, da seguinte forma. Como um exemplo, nós ilustramos uma
tentativa de adicionar ao preambulo um pacote que supostamente não
existe. ::

    sage: latex.extra_preamble('')
    sage: latex.extra_preamble()
    ''
    sage: latex.add_to_preamble('\\usepackage{foo-bar-unchecked}')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'
    sage: latex.add_package_to_preamble_if_available('foo-bar-checked')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'

.. _sec-custom-processing:

Personalizando o Processamento em LaTeX
=======================================

É também possível controlar qual variação do TeX é usada quando a
versão do sistema for evocada, logo influenciando também o resultado.
De forma similar, é também possível controlar quando o Notebook irá
usar o MathJax (trechos simples em TeX) ou a versão do TeX do sistema
(expressões mais complicadas).

O comando ``latex.engine()`` pode ser usado para controlar de os
executáveis ``latex``, ``pdflatex`` ou ``xelatex`` do sistema são
usados para processar expressões mais complicadas. Quando ``view()`` é
chamado na linha de comando do Sage e o processador é definido como
``latex``, um arquivo dvi é produzido e o Sage vai usar um
visualizador de dvi (como o xdvi) para apresentar o resultado. Por
outro lado, usando ``view()`` na linha de comando do Sage, quando o
processador é definido como ``pdflatex``, irá produzir um PDF e o Sage vai
executar o programa disponível no seu sistema para visualizar arquivos
PDF (acrobat, okular, evince, etc.).

No Notebook, é necessário interver na decisão de se o MathJax vai
interpretar trechos em TeX, ou se o LaTeX do sistema deve fazer o
trabalho se o código em LaTeX for complicado demais. O dispositivo é
uma lista de strings, que se forem encontradas em um trecho de código
LaTeX sinalizam para o Notebook usar o LaTeX (ou qualquer executável
que for definido pelo comando ``latex.engine()``). Essa lista é
gerenciada pelos comandos ``latex.add_to_mathjax_avoid_list`` e
``latex.mathjax_avoid_list``. ::

    sage: # not tested
    sage: latex.mathjax_avoid_list([])
    sage: latex.mathjax_avoid_list()
    []
    sage: latex.mathjax_avoid_list(['foo', 'bar'])
    sage: latex.mathjax_avoid_list()
    ['foo', 'bar']
    sage: latex.add_to_mathjax_avoid_list('tikzpicture')
    sage: latex.mathjax_avoid_list()
    ['foo', 'bar', 'tikzpicture']
    sage: latex.mathjax_avoid_list([])
    sage: latex.mathjax_avoid_list()
    []

Suponha que uma expressão em LaTeX é produzida no Notebook com o
comando ``view()`` ou enquanto o botão "Typeset" está selecionado, e
então reconhecida, através da "lista de comandos a serem evitados no
MathJax", como necessitando a versão do LaTeX no sistema. Então o
executável selecionado (como especificado por ``latex.engine()``) irá
processar o código em LaTeX. Todavia, em vez de então abrir um
visualizador externo (o que é o comportamento na linha de comando), o
Sage irá tentar converter o resultado em uma imagem, que então é
inserida na folha de trabalho como o resultado da célula.

Exatamente como essa conversão é feita depende de vários fatores --
qual executável você especificou como processador e quais utilitários
de conversão estão disponíveis no seu sistema. Quatro conversores
usuais que irão cobrir todas as ocorrências são o ``dvips``,
``ps2pdf``, e ``dvipng``, e do pacote ``ImageMagick``, o ``convert``.
O objetivo é produzir um arquivo PNG para ser inserido de volta na
folha de trabalho. Quando uma expressão em LaTeX pode ser convertida
com sucesso em um arquivo dvi pelo processador LaTeX, então o dvipng
deve dar conta da conversão. Se a expressão em LaTeX e o processador
especificado criarem um arquivo dvi com conteúdo especial que o dvipng
não pode converter, então o dvips vai criar um arquivo PostScript.
Esse arquivo PostScript, ou um PDF criado por pelo processador
``pdflatex``, é então convertido em um arquivo dvi pelo programa
``convert``. A presença de dois desses conversores pode ser testado
com as rotinas ``have_dvipng()`` e ``have_convert()``.

Essas conversões são feitas automaticamente se você tiver os
conversores necessários instalados; se não, então uma mensagem de erro
é impressa dizendo o que está faltando e onde obter.

Para um exemplo concreto de como expressões complicadas em LaTeX podem
ser processadas, veja o exemplo na próxima seção
(:ref:`sec-tkz-graph`) para usar o pacote ``tkz-graph`` para produzir
ilustrações de grafos combinatoriais de alta qualidade. Para outros
exemplos, existem alguns casos teste incluídos no Sage. Para usá-los,
é necessário importar o objeto ``sage.misc.latex.latex_examples``, que
é uma instância da classe ``sage.misc.latex.LatexExamples``, como
mostrado abaixo. Essa classe possui exemplos de diagramas comutativos,
grafos combinatoriais, teoria de nós e pstricks, os quais
respectivamente testam os seguintes pacotes: xy, tkz-graph, xypic,
pstricks. Após importar o objeto, use completamento tab em
``latex_examples`` para ver os exemplos disponíveis. Ao carregar um
exemplo você irá obter explicações sobre o que é necessário para fazer
o conteúdo do exemplo ser exibido corretamente. Para de fato ver os
exemplos, é necessário usar ``view()`` (uma vez que o preambulo,
processador, etc. estão configurados corretamente).

::

    sage: from sage.misc.latex import latex_examples
    sage: latex_examples.diagram()
    LaTeX example for testing display of a commutative diagram produced
    by xypic.
    <BLANKLINE>
    To use, try to view this object -- it will not work.  Now try
    'latex.add_to_preamble("\\usepackage[matrix,arrow,curve,cmtip]{xy}")',
    and try viewing again. You should get a picture (a part of the diagram arising
    from a filtered chain complex).

.. _sec-tkz-graph:

Exemplo: Grafos Combinatoriais com tkz-graph
============================================

Ilustrações de alta qualidade de grafos combinatoriais (daqui por
diante, simplesmente grafos) são possíveis com o pacote ``tkz-graph``.
Esse pacote baseia-se no ``tikz`` front-end da biblioteca ``pgf``.
Logo todos esses componentes precisam ser parte de uma instalação
completa do LaTeX em seu sistema, e pode acontecer que alguns desses
componentes não estejam em sua versão mais recente em algumas
distribuições do TeX. Logo, para melhores resultados, seria necessário
ou recomendável instalar esses pacotes como parte do seu diretório
texmf pessoal. Criar, manter e personalizar uma instalação do TeX no
sistema ou em um diretório pessoal vai além do escopo deste documento,
mas deve ser fácil encontrar instruções para isso. Os arquivos
necessários estão listados em :ref:`sec-system-wide-tex`.

Portanto, para começar precisamos nos certificar que os pacotes
relevantes estão incluídos adicionando-os ao preambulo do eventual
documento LaTeX. As imagens dos grafos não são formadas corretamente
quando um arquivo dvi é usando como formato intermediário, logo é
melhor definir o processador do LaTeX como ``pdflatex``. A esta altura
um comando como ``view(graphs.CompleteGraph(4))`` deve funcionar na
linha de comando do Sage e produzir um PDF com a imagem completa do
grafo `K_4`.

Para uma experiência semelhante no Notebook, é necessário desabilitar
o processador MathJax para o código LaTeX do grafo usando a "lista de
comandos a serem evitados pelo MathJax". Grafos são criados usando o
ambiente ``tikzpicture``, logo essa uma boa escolha para uma string
a ser incluída na lista que acabamos de mencionar. Agora,
``view(graphs.CompleteGraph(4))`` em uma folha de trabalho deve
executar o pdflatex para criar um PDF e então o programa ``convert``
para obter um gráfico PNG que vai ser inserido na folha de trabalho.
Os seguintes comandos ilustram os passos para obter grafos processados
pelo LaTeX no Notebook. ::

    sage: from sage.graphs.graph_latex import setup_latex_preamble
    sage: setup_latex_preamble()
    sage: latex.extra_preamble() # random - depends on system's TeX installation
    '\\usepackage{tikz}\n\\usepackage{tkz-graph}\n\\usepackage{tkz-berge}\n'
    sage: latex.engine('pdflatex')
    sage: latex.add_to_mathjax_avoid_list('tikzpicture')  # not tested
    sage: latex.mathjax_avoid_list()                      # not tested
    ['tikz', 'tikzpicture']

Agora, um comando como ``view(graphs.CompleteGraph(4))`` deve produzir
um gráfico do grafo no Notebook, tendo usado ``pdflatex`` para
processar os comandos do ``tkz-graph`` para construir o grafo. Note
que há diversas opções que afetam o resultado do gráfico obtido usando
o LaTeX via ``tkz-graph``, o que mais uma vez está além do escopo
desta seção (veja a seção do Manual de Referência com título "Opções
do LaTeX para Grafos" para instruções e detalhes).

.. _sec-system-wide-tex:

Uma Instalação Completa do TeX
==============================
Vários dos recursos avançados de integração do TeX com o Sage requerem
uma instalação do TeX em seu sistema. Várias versões do Linux possuem
pacotes do TeX baseados no TeX-live, para o OSX existe o TeXshop e
para o windows existe o MikTex. O utilitário ``convert`` é parte do
`ImageMagick <http://www.imagemagick.org/>`_ (que deve ser um pacote
na sua versão do Linux ou ser fácil de instalar), e os três programas
``dvipng``, ``ps2pdf``, e ``dvips`` podem estar incluídos na sua
distribuição do TeX. Os dois primeiros podem também ser obtidos em,
respectivamente, http://sourceforge.net/projects/dvipng/ e como parte
do `Ghostscript <http://www.ghostscript.com/>`_.

A criação de grafos combinatoriais requer uma versão recente da
biblioteca PGF, e o arquivo ``tkz-graph.sty``, que está disponível em
https://www.ctan.org/pkg/tkz-graph, os arquivos  ``tkz-arith.sty`` e
talvez ``tkz-berge.sty``, que estão disponíveis em
https://www.ctan.org/pkg/tkz-berge.

Programas Externos
==================

Existem três programas disponíveis para integrar ainda mais o TeX e o
Sage. O primeiro é o sagetex. Uma descrição concisa do sagetex é que
ele é uma coleção de funções do TeX que permitem incluir em um
documento LaTeX instruções para usar o Sage para calcular vários
objetos, e/ou formatar objetos usando o comando ``latex()`` existente
no Sage. Logo, como um passo intermediário para compilar um documento
LaTeX, todos os recursos computacionais e de formatação do Sage podem
ser executados automaticamente. Como um exemplo, um exame matemático
pode manter uma correspondência entre questões e respostas usando o
sagetex para fazer cálculos com o Sage. Veja :ref:`sec-sagetex` para
mais informações.
