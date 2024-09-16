.. _chapter-interactive_shell:

*****************************
A Linha de Comando Interativa
*****************************
Na maior parte deste tutorial, assumimos que você iniciou o
interpretador Sage usando o comando ``sage``. Isso inicia uma versão
personalizada da linha de comando IPython, e importa diversas funções
e classes de modo que elas fiquem prontas para serem usadas a partir
da linha de comando. Configuração adicional é possível editando o
arquivo ``$SAGE_ROOT/ipythonrc``. Assim que você inicia o Sage, você
obtém o seguinte:

.. skip

::

    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘
    sage:

Para sair do Sage pressione Ctrl-D ou digite ``quit`` ou ``exit``.

.. skip

::

    sage: quit
    Exiting Sage (CPU time 0m0.00s, Wall time 0m0.89s)

O wall time é o tempo que passou no relógio "pendurado na sua parede".
Isso é relevante, pois o tempo CPU não conta o tempo usado por
subprocessos como GAP ou Singular.

(Evite terminar um processo do Sage usando ``kill -9`` a partir de um
terminal, pois o Sage pode não terminal seus subprocessos, por
exemplo, subprocessos do Maple, ou limpeza de arquivos temporários em
``$HOME/.sage/tmp``.)

A Sua Sessão no Sage
====================

A sessão é a sequência de entradas e saídas de dados desde o momento em
que você inicia até o momento em que você termina o Sage. O Sage grava
todas as entradas de dados, através do IPython. De fato, se você está
usando a linha de comando (não o Notebook), então a qualquer momento
você pode digitar ``%history`` (ou ``%hist``) para obter uma lista de
todas as linhas digitadas até então. Você pode digitar ``?`` no prompt
do Sage para aprender mais sobre o IPython, por exemplo, "IPython
offers numbered prompts ... with input and output caching. All input
is saved and can be retrieved as variables (besides the usual arrow
key recall). The following GLOBAL variables always exist (so don't
overwrite them!)":

::

      _:  previous input (interactive shell and notebook)
      __: next previous input (interactive shell only)
      _oh : list of all inputs (interactive shell only)

Aqui vai um exemplo:

.. skip

::

    sage: factor(100)
     _1 = 2^2 * 5^2
    sage: kronecker_symbol(3,5)
     _2 = -1
    sage: %hist   #This only works from the interactive shell, not the notebook.
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    sage: _oh
     _4 = {1: 2^2 * 5^2, 2: -1}
    sage: _i1
     _5 = 'factor(ZZ(100))\n'
    sage: eval(_i1)
     _6 = 2^2 * 5^2
    sage: %hist
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    4: _oh
    5: _i1
    6: eval(_i1)
    7: %hist

Vamos omitir a numeração das linhas no restante deste tutorial e em
outras documentações do Sage.

Você também pode salvar uma lista de comandos em uma macro.

.. skip

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: M = ModularSymbols(37)
    sage: %hist
    1: E = EllipticCurve([1,2,3,4,5])
    2: M = ModularSymbols(37)
    3: %hist
    sage: %macro em 1-2
    Macro `em` created. To execute, type its name (without quotes).


.. skip

::

    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E = 5
    sage: M = None
    sage: em
    Executing Macro...
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field

Quando se usa a linha de comando, qualquer comando UNIX pode ser
executado a partir do Sage inserindo um ponto de exclamação ``!`` como
prefixo. Por exemplo,

.. skip

::

    sage: !ls
    auto  example.sage glossary.tex  t  tmp  tut.log  tut.tex

fornece a lista de arquivos do atual diretório.

O ``PATH`` possui o diretório bin do Sage em primeiro, portanto se
você digitar ``p``, ``gap``, ``singular``, ``maxima``, etc., você
executa a versão incluída no Sage.

.. skip

::

    sage: !gp
    Reading GPRC: /etc/gprc ...Done.

                               GP/PARI CALCULATOR Version 2.2.11 (alpha)
                      i686 running linux (ix86/GMP-4.1.4 kernel) 32-bit version
    ...
    sage: !singular
                         SINGULAR                             /  Development
     A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \

Gravando Entradas e Saídas de dados
===================================

Gravar a sua sessão no Sage não é o mesmo que salvá-la (veja
:ref:`section-save`). Para gravar a entrada de dados (e opcionalmente
a saída) use o comando ``logstart``. Digite ``logstart?`` para mais
detalhes. Você pode usar esse comando para gravar tudo o que você
digita, toda a saída de dados, e até mesmo usar essa entrada de dados
que você guardou em uma sessão futura (simplesmente importando o
arquivo log).

.. skip

::

    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: logstart setup
    Activating auto-logging. Current session state plus future input saved.
    Filename       : setup
    Mode           : backup
    Output logging : False
    Timestamping   : False
    State          : active
    sage: E = EllipticCurve([1,2,3,4,5]).minimal_model()
    sage: F = QQ^3
    sage: x,y = QQ['x,y'].gens()
    sage: G = E.gens()
    sage:
    Exiting Sage (CPU time 0m0.61s, Wall time 0m50.39s).
    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: load "setup"
    Loading log file <setup> one line at a time...
    Finished replaying log file <setup>
    sage: E
    Elliptic Curve defined by y^2 + x*y  = x^3 - x^2 + 4*x + 3 over Rational
    Field
    sage: x*y
    x*y
    sage: G
    [(2 : 3 : 1)]

Se você usa o Sage no terminal ``konsole`` do Linux KDE, então você
pode gravar a sessão da seguinte forma: após iniciar o Sage no
``konsole``, selecione "settings", então "history...", então "set
unlimited". Quando você estiver pronto para guardar a sua sessão,
selecione "edit" e então "save history as..." e digite um nome para
salvar o texto de sua sessão em seu computador. Após salvar esse
arquivo, você poderia abri-lô em um editor, tal como o xemacs, e
imprimi-lo.

Colar Texto Ignora Prompts
==========================

Suponha que você está lendo uma sequência de comandos em Sage ou
Python e quer copiá-los no Sage. Mas eles têm os irritantes prompts
``>>>`` ou ``sage:`` para te aborrecer. De fato, você pode copiar e
colar um exemplo, incluindo os prompts se você quiser, no Sage. Em
outras palavras, automaticamente o Sage remove os caracteres ``>>>``
ou ``sage:`` antes de colar o conteúdo no Python. Por exemplo,

.. skip

::

    sage: 2^10
    1024
    sage: sage: sage: 2^10
    1024
    sage: >>> 2^10
    1024

Comandos de Tempo
=================

Se você colocar o comando ``%time`` no começo de uma linha de comando,
o tempo que o comando leva para ser executado vai aparecer após a
saída de dados. Por exemplo, nós podemos comparar o tempo de execução
para certas operações de exponenciação de várias formas. Os tempos
abaixo vão ser provavelmente muito diferentes para o seu computador,
ou até mesmo para versões diferentes do Sage. Primeiro, usando o
Python

.. skip

::

    sage: %time a = int(1938)^int(99484)
    CPU times: user 0.66 s, sys: 0.00 s, total: 0.66 s
    Wall time: 0.66

Isso significa que levou 0.66 segundos no total, e o wall time, isto
é, a quantidade de tempo que passou no seu relógio de parede, é também
0.66 segundos. Se o seu computador está executado outros programas o
wall time pode ser muito maior do que o tempo de CPU.

A seguir verificamos o tempo de exponenciação usando o tipo Integer do
Sage, o qual é implementado (em Cython) usando a biblioteca GMP:

.. skip

::

    sage: %time a = 1938^99484
    CPU times: user 0.04 s, sys: 0.00 s, total: 0.04 s
    Wall time: 0.04

Usando a biblioteca C do PARI:

.. skip

::

    sage: %time a = pari(1938)^pari(99484)
    CPU times: user 0.05 s, sys: 0.00 s, total: 0.05 s
    Wall time: 0.05

A GMP é melhor, mas por pouco (como esperado, pois a versão do PARI
contida no Sage usa a GMP para aritmética de inteiros).

Você pode também contar o tempo de um bloco de comandos usado o
comando ``cputime``, como ilustrado abaixo:

::

    sage: t = cputime()
    sage: a = int(1938)^int(99484)
    sage: b = 1938^99484
    sage: c = pari(1938)^pari(99484)
    sage: cputime(t)                       # somewhat random output
    0.64

.. skip

::

    sage: cputime?
    ...
        Return the time in CPU second since Sage started, or with optional
        argument t, return the time since time t.
        INPUT:
            t -- (optional) float, time in CPU seconds
        OUTPUT:
            float -- time in CPU seconds

O comando ``walltime`` se comporta como o comando ``cputime``, exceto
que ele conta o tempo do relógio.

Nós podemos também calcular a potência acima em alguns softwares de
álgebra incluídos no Sage. Em cada caso executamos um comando trivial
no sistema de modo a inicializar o servidor para aquele programa. O
tempo mais relevante é o tempo do relógio. Todavia, se houver uma
diferença significativa entre o wall time e o CPU time então isso pode
indicar algum problema de performance que vale a pena investigar.

.. skip

::

    sage: time 1938^99484;
    CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
    Wall time: 0.01
    sage: gp(0)
    0
    sage: time g = gp('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: maxima(0)
    0
    sage: time g = maxima('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.30
    sage: kash(0)
    0
    sage: time g = kash('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: mathematica(0)
            0
    sage: time g = mathematica('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.03
    sage: maple(0)
    0
    sage: time g = maple('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.11
    sage: libgap(0)
    0
    sage: time g = libgap.eval('1938^99484;')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 1.02

Note que o GAP e o Maxima são os mais lentos neste teste (isso foi
executado no computador ``sage.math.washington.edu``). Devido ao
processamento extra (overhead) da interface pexpect, talvez não seja
apropriado comparar esses resultados com o Sage, que é o mais rápido.

Outras Dicas para o IPython
===========================

Como observado acima, o Sage usa o IPython como interface, logo você
pode usar quaisquer comandos e recursos do IPython. Você pode ler a
`Documentação completa do IPython
<http://ipython.scipy.org/moin/Documentation>`_ (em inglês).

- Você pode usar ``%edit`` (ou ``%ed`` ou ``ed``) para abrir um
  editor, se você desejar digitar algum código mais complexo. Antes de
  iniciar o Sage, certifique-se de que a variável de ambiente
  :envvar:`EDITOR` está definida com o seu editor favorito (colocando
  ``export EDITOR=/usr/bin/emacs`` ou ``export EDITOR=/usr/bin/vim``
  or algo similar no lugar apropriado, como um arquivo ``.profile``).
  A partir do prompt do Sage, o comando ``%edit`` irá abrir o editor
  escolhido. Então usando o editor você pode definir uma função::

    def some_function(n):
        return n**2 + 3*n + 2

  Salve e feche o editor. No restante da sua sessão do Sage, você pode
  usar então a função ``some_function``. Se você quiser modificá-la,
  digite ``%edit some_function`` no prompt do Sage.

- Se você for calcular algo e quiser modificar o resultado para outro
  uso, execute o cálculo e então digite ``%rep``: isso irá colocar o
  resultado do comando anterior no prompt do Sage, pronto para ser
  editado.::

    sage: f(x) = cos(x)
    sage: f(x).derivative(x)
    -sin(x)

  A esta altura, se você digitar ``%rep`` no prompt do Sage, você irá
  obter um novo prompt, seguido de ``-sin(x)``, com o cursor no final
  da linha.

Para mais informações, digite ``%quickref`` para ver um guia rápido de
referência do IPython. Quando este tutorial foi escrito (Fevereiro 2015),
o Sage usa a versão 2.3.0 do IPython, e a `documentation for its magic commands <http://ipython.org/ipython-doc/dev/interactive/tutorial.html#magic-functions>`_ está disponível na internet.


Erros e Exceções
================

Quando algo errado ocorre, você usualmente verá uma "exceção" do
Python. O Python até mesmo tenta sugerir o que ocasionou a exceção,
por exemplo, :class:`NameError` ou :class:`ValueError` (veja o Referência da
Biblioteca Python [PyLR]_ para uma lista completa de exceções). Por
exemplo,

.. skip

::

    sage: 3_2
    ------------------------------------------------------------
       File "<console>", line 1
         ZZ(3)_2
               ^
    SyntaxError: invalid ...

    sage: EllipticCurve([0,infinity])
    ------------------------------------------------------------
    Traceback (most recent call last):
    ...
    TypeError: Unable to coerce Infinity (<class 'sage...Infinity'>) to Rational

O debugador interativo é as vezes útil para entender o que houve de
errado. Você pode ativá-lo e desativá-lo usando ``%pdb`` (o padrão é
desativado). O prompt ``ipdb>`` aparece se uma exceção é levantada e o
debugador está ativado. A partir do debugador, você pode imprimir o
estado de qualquer variável local, e mover a pilha de execução para
cima e para baixo. Por exemplo,

.. skip

::

    sage: %pdb
    Automatic pdb calling has been turned ON
    sage: EllipticCurve([1,infinity])
    ---------------------------------------------------------------------------
    <class 'exceptions.TypeError'>             Traceback (most recent call last)
    ...

    ipdb>

Para uma lista de comandos do debugador, digite ``?`` no prompt
``ipbd>``:

::

    ipdb> ?

    Documented commands (type help <topic>):
    ========================================
    EOF    break  commands   debug    h       l     pdef   quit    tbreak
    a      bt     condition  disable  help    list  pdoc   r       u
    alias  c      cont       down     ignore  n     pinfo  return  unalias
    args   cl     continue   enable   j       next  pp     s       up
    b      clear  d          exit     jump    p     q      step    w
    whatis where

    Miscellaneous help topics:
    ==========================
    exec  pdb

    Undocumented commands:
    ======================
    retval  rv

Digite Ctrl-D ou ``quit`` para retornar ao Sage.

.. _section-tabcompletion:

Busca Reversa e Completamento Tab
==================================

Busca reversa: Digite o começo de um comando, e então ``Ctrl-p`` (ou
tecle a seta para cima) para voltar para cada linha que você digitou
que começa daquela forma. Isso funciona mesmo se você encerrou o Sage
e iniciou novamente mais tarde. Você também pode fazer uma busca
reversa ao longo da história usando ``Ctrl-r``. Todos esses recursos
usam o pacote ``readline``, que está disponível no Linux.

Para ilustrar a busca reversa, primeiro crie o e espaço vetorial
tri-dimensional :math:`V=\QQ^3` da seguinte forma:

::

    sage: V = VectorSpace(QQ,3)
    sage: V
    Vector space of dimension 3 over Rational Field

Você pode usar a seguinte notação mais compacta:

::

    sage: V = QQ^3

Então é fácil listar todas as funções para :math:`V` usando
completamento. Digite ``V``, e então pressione a tecla ``[tab]`` no
seu teclado:

.. skip

::

    sage: V.[tab key]
    V._VectorSpace_generic__base_field
    ...
    V.ambient_space
    V.base_field
    V.base_ring
    V.basis
    V.coordinates
    ...
    V.zero_vector

Se você digitar as primeiras letras de uma função, e então a tecla
``[tab]``, você obtém apenas funções que começam conforme indicado.

.. skip

::

    sage: V.i[tab key]
    V.is_ambient  V.is_dense    V.is_full     V.is_sparse

Se você gostaria de saber o que uma função faz, por exemplo, a função
coordinates, digite ``V.coordinates?`` para ajuda ou
``V.coordinates??`` para ver o código fonte, como explicado na próxima
sessão.


Sistema de Ajuda Integrado
==========================

O Sage possui um sistema de ajuda integrado. Digite o nome da função
seguido de ? para ver informações sobre a função.

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates?
    Type:           instancemethod
    Base Class:     <class 'instancemethod'>
    String Form:    <bound method FreeModule_ambient_field.coordinates of Vector
    space of dimension 3 over Rational Field>
    Namespace:      Interactive
    File:           /home/was/s/local/lib/python2.4/site-packages/sage/modules/f
    ree_module.py
    Definition:     V.coordinates(self, v)
    Docstring:
        Write v in terms of the basis for self.

        Returns a list c such that if B is the basis for self, then

                sum c_i B_i = v.

        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: M = FreeModule(IntegerRing(), 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0-M1)
            [2, -1]

Como mostrado acima, o comando de ajuda mostra o tipo do objeto, o
arquivo no qual ele é definido, e uma descrição útil da função com
exemplos que você pode colar na sua sessão atual. Quase todos esses
exemplos são regularmente testados automaticamente para certificar que
eles se comportam exatamente como esperado.

Outro recurso que vai muito na direção do espírito de código aberto do
Sage é que se ``f`` é uma função do Python, então o comando ``f??``
mostra o código fonte que define ``f``. Por exemplo,

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates??
    Type:           instancemethod
    ...
    Source:
    def coordinates(self, v):
            """
            Write $v$ in terms of the basis for self.
            ...
            """
            return self.coordinate_vector(v).list()

Isso nos diz que tudo que a função ``coordinates`` faz é chamar a
função ``coordinate_vector`` e converter o resultado para uma lista. O
que a função ``coordinate_vector`` faz?

.. skip

::

    sage: V = QQ^3
    sage: V.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            ...
            return self.ambient_vector_space()(v)

A função ``coordinate_vector`` coage a sua entrada em um espaço
ambiente, o que tem o efeito de calcular o vetor de coeficientes de
:math:`v` em termos de :math:`V`. O espaço :math:`V` já é o espaço
ambiente pois é simplesmente :math:`\QQ^3`. Existe também uma função
``coordinate_vector`` para subespaços, que é diferente. Vamos criar um
subespaço e ver:

.. skip

::

    sage: V = QQ^3; W = V.span_of_basis([V.0, V.1])
    sage: W.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            """
             ...
            """
            # First find the coordinates of v wrt echelon basis.
            w = self.echelon_coordinate_vector(v)
            # Next use transformation matrix from echelon basis to
            # user basis.
            T = self.echelon_to_user_matrix()
            return T.linear_combination_of_rows(w)

(Se você acha que a implementação é ineficiente, por favor junte-se a
nós para ajudar a optimizar as funções de álgebra linear.)

Você também pode digitar ``help(command_name)`` ou ``help(class)``
para ler um arquivo de ajuda sobre determinada classe.

.. skip

::

    sage: help(VectorSpace)
    Help on function VectorSpace in module sage.modules.free_module:

    VectorSpace(K, dimension_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
                with_basis='standard', dimension=None, basis_keys=None, **args)
    EXAMPLES:

    The base can be complicated, as long as it is a field.

    ::

        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x
         over Integer Ring
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
    --More--

Quando você digita ``q`` para sair do sistema de ajuda, a sua sessão
aparece na tela da mesma forma que anteriormente. O texto de ajuda não
fica permanentemente em sua tela, ao contrário da saída de
``function_name?`` que as vezes fica. É partircularmente útil digitar
``help(module_name)``. Por exemplo, espaços vetoriais são definidos em
``sage.modules.free_module``, então digite
``help(sage.modules.free_module)`` para obter documentação sobre esse
módulo. Quando visualizando documentação usando a ajuda, você pode
procurar no texto digitando ``/`` e na ordem reversa digitando ``?``.

Salvando e Carregando Objetos Individuais
=========================================

Suponha que você calcule uma matriz, ou pior ainda, um espaço
complicado de símbolos, e gostaria de salvá-los para uso posterior. O
que você pode fazer? Existem várias estratégias que os sistemas
computacionais de álgebra adotam para salvar objetos individuais.

#. **Salve seus cálculos:** Suportam apenas salvar e carregar uma
   sessão completa (por exemplo, GAP, Magma).

#. **Entrada e Saída Unificadas:** Faz com que cada objeto seja
   impresso de uma forma que possa ser lido novamente (GP/PARI).

#. **Eval:** Torna fácil processar um código arbitrário no
   interpretador (por exemplo, Singular, PARI).


Como o Sage usa o Python, ele adota uma estratégia diferente, que se
baseia no fato de que qualquer objeto pode ser "serializado", isto é,
transformado em uma string a partir da qual o objeto pode ser
recuperado. Isso segue o espírito da estratégia unificada de entrada e
saída do PARI, exceto que não possue a desvantagem que os objetos são
impressos na tela em uma forma muito complicada. Além disso, o suporte
para salvar e recuperar é (na maior parte dos casos) completamente
automática, não requerendo nenhuma programação extra; é simplesmente um
recurso do Python que foi implementado na linguagem desde o início de
seu desenvolvimento.

Quase todos os objetos ``x`` podem ser armazenados em disco de forma
comprimida usando ``save(x, filename)`` (ou em muitos casos
``x.save(filename)``). Para carregar o objeto de volta no Sage use
``load(filename)``.

.. skip

::

    sage: A = MatrixSpace(QQ,3)(range(9))^2
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]
    sage: save(A, 'A')

Você deve agora sair do Sage e reiniciá-lo. Então você pode obter
``A`` de volta:

.. skip

::

    sage: A = load('A')
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]

Você pode fazer o mesmo com objetos mais complicados, por exemplo,
curvas elípticas. Todos os dados sobre o objeto são guardados e
restaurados com o objeto. Por exemplo,

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: v = E.anlist(100000)              # takes a while
    sage: save(E, 'E')
    sage: quit

A versão em disco de ``E`` ocupa 153 kilobytes, pois ela guarda os
primeiros 1000000 :math:`a_n` com ela.

.. skip

::

    ~/tmp$ ls -l E.sobj
    -rw-r--r--  1 was was 153500 2006-01-28 19:23 E.sobj
    ~/tmp$ sage [...]
    sage: E = load('E')
    sage: v = E.anlist(100000)              # instant!

(Em Python, salvar e restaurar é feito usando o módulo ``cPickle``. Em
particular, um objeto ``x`` do Sage pode ser salvo usando
``cPickle.dumps(x, 2)``. Note o ``2``!)

O sage não pode salvar e carregar objetos criados em algum outro
sistema computacional de álgebra, por exemplo, GAP, Singular, Maxima,
etc. Eles são carregados em um estado "inválido". Em GAP, embora
muitos objetos podem ser impressos de uma forma que eles podem ser
reconstruídos, muitos não, logo reconstrução a partir de suas
representações impressas não é permitido.

.. skip

::

    sage: a = libgap(2)
    sage: a.save('a')
    sage: load('a')
    Traceback (most recent call last):
    ...
    ValueError: The session in which this object was defined is no longer
    running.

Objetos do GP/PARI também podem ser salvos e carregados pois suas
representações em forma impressa são suficientes para reconstruí-los.

.. skip

::

    sage: a = gp(2)
    sage: a.save('a')
    sage: load('a')
    2

Objetos que foram salvos podem ser abertos posteriormente em
computadores com arquiteturas e sistemas operacionais diferentes, por
exemplo, você poderia salvar uma matriz muito grande em um OS X de
32-bits e abri-lo em um Linux de 64-bits, encontrar a forma reduzida,
e movê-lo de volta. Além disso, em muitos casos você pode até mesmo
abrir objetos em versões do Sage diferentes daquela no qual o objeto
foi salvo, desde que o código para aquele objeto não seja muito
diferente. Todos os atributos do objetos são armazendos, assim como a
classe (mas não o código fonte) que define o objeto. Se aquela classe
não existir mais em uma nova versão do Sage, então o objeto não pode
ser reaberto nessa versão. Mas você poderia abri-lo em uma versão mais
antiga, obter o dicionário do objeto (com ``x.__dict__``), salvar o
dicionário, e abri-lo em uma versão mais nova.

Salvando como Texto
-------------------

Você também pode salvar a representação em texto (ASCII) de objetos em
um arquivo, simplesmente abrindo um arquivo em modo de escrita, e
escrevendo a string que representa o objeto no arquivo (você pode
salvar mais de um objeto dessa forma). Quando você terminar de
escrever os objetos, feche o arquivo.

.. skip

::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: f = (x+y)^7
    sage: o = open('file.txt','w')
    sage: o.write(str(f))
    sage: o.close()

.. _section-save:

Salvando e Abrindo Sessões Completas
====================================

O Sage é flexível para salvar e abrir sessões completas.

O comando ``save_session(sessionname)`` salva todas as variáveis que
você definiu na sessão atual como um dicionário com o nome
``sessionname``. (No caso raro de uma variável não poder ser salva,
ela simplesmente não aparece no dicionário.) O resultado é um arquivo
``.sobj`` que pode ser aberto como qualquer outro objeto que foi
salvo. Quando você abre os objetos que foram salvos em uma sessão,
você obtém um dicionário cujas chaves (keys) são os nomes das
variáveis e os valores são os objetos.

Você pode usar o comando ``load_session(sessionname)`` para carregar
na presente sessão as variáveis definidas em ``sessioname``. Note que
isso não remove as variáveis já definidas na presente sessão; em vez
disso, as duas sessões são combinadas.

Primeiro iniciamos o Sage e definimos algumas variáveis.

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: M = ModularSymbols(37)
    sage: a = 389
    sage: t = M.T(2003).matrix(); t.charpoly().factor()
     _4 = (x - 2004) * (x - 12)^2 * (x + 54)^2

A seguir nós salvamos a nossa sessão, o que armazena cada uma das
variáveis acima em um arquivo. Então visualizamos o arquivo, que tem
por volta de 3K bytes.

.. skip

::

    sage: save_session('misc')
    Saving a
    Saving M
    Saving t
    Saving E
    sage: quit
    was@form:~/tmp$ ls -l misc.sobj
    -rw-r--r--  1 was was 2979 2006-01-28 19:47 misc.sobj

Por fim reiniciamos o Sage, definimos algumas variáveis extra, e
carregamos a sessão que foi salva anteriormente.

.. skip

::

    sage: b = 19
    sage: load_session('misc')
    Loading a
    Loading M
    Loading E
    Loading t

Cada variável que foi salva está de novo disponível. Além disso, a
variável ``b`` não foi redefinida.

.. skip

::

    sage: M
    Full Modular Symbols space for Gamma_0(37) of weight 2 with sign 0
    and dimension 5 over Rational Field
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational
    Field
    sage: b
    19
    sage: a
    389

