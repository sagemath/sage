***********
Programação
***********

.. _section-loadattach:

Carregando e Anexando Arquivos do Sage
======================================

A seguir ilustramos como carregar no Sage programas escritos em um
arquivo separado. Crie um arquivo chamado ``example.sage`` com o
seguinte conteúdo:

.. skip

::

    print("Hello World")
    print(2^3)

Você pode ler e executar o arquivo ``example.sage`` usando o comando
``load``.

.. skip

::

    sage: load "example.sage"
    Hello World
    8

Você também pode anexar um arquivo em Sage à sessão em execução usando
o comando ``attach``.

.. skip

::

    sage: attach "example.sage"
    Hello World
    8

Agora se você alterar ``example.sage`` e adicionar uma linha em branco
(por exemplo), então o conteúdo de ``example.sage`` será
automaticamente recarregado no Sage.

Em particular, ``attach`` automaticamente recarrega um arquivo toda
vez que ele for modificado, o que é útil para desenvolver e testar um
programa, enquanto ``load`` carrega o arquivo apenas uma vez.

Quando o Sage carrega ``example.sage`` ele converte esse arquivo para
o Python, o qual é então executado pelo interpretador do Python. Essa
conversão é mínima; ela essencialmente consiste em encapsular inteiros
em ``Integer()``, números float em ``RealNumber()``, substituir ``^``
por ``**``, e substituir, por exemplo, ``R.2`` por ``R.gen(2)``. A
versão convertida de ``example.sage`` é armazenada no mesmo diretório
de ``example.sage`` e é chamada ``example.sage.py``. Esse arquivo
contém o seguinte código:

::

    print("Hello World")
    print(Integer(2)**Integer(3))

Inteiros literais são encapsulados e ``^`` é substituído por ``**``.
(Em Python, ``^`` significa "ou exclusivo" e ``**`` significa
"exponenciação".)

Esse "" está implementado em ``sage/misc/interpreter.py``.)

Você pode colar código tabulado com muitas linhas no Sage desde que
existam linhas em branco separando blocos de código (isso não é
necessário em arquivos). Todavia, a melhor forma de adicionar tal
código a uma sessão do Sage é salvá-lo em um arquivo e usar
``attach``, como descrito anteriormente.


.. _section-compile:

Criando Código Compilado
========================

Velocidade é crucial em cálculos matemáticos. Embora o Python seja uma
linguagem conveniente de alto nível, certos cálculos podem ser várias
vezes mais rápidos do que em Python se eles forem implementados usando
tipos estáticos em uma linguagem compilada. Alguns aspectos do Sage
seriam muito lentos se eles fossem escritos inteiramente em Python.
Para lidar com isso, o Sage suporta uma "versão" compilada do Python
chamada Cython ([Cyt]_ and [Pyr]_). O Cython é simultaneamente similar
ao Python e ao C. A maior parte das construções em Python, incluindo
"list comprehensions", expressões condicionais, código como ``+=`` são
permitidos; você também pode importar código que você escreveu em
outros módulos em Python. Além disso, você pode declarar variáveis em
C arbitrárias, e qualquer chamada de bibliotecas em C pode ser feita
diretamente. O código resultante é convertido para C e compilado
usando um compilador para C.

Para criar o seu próprio código compilado em Sage, nomeie o arquivo
com uma extensão ``.spyx`` (em vez de ``.sage``). Se você está
trabalhando na linha de comando, você pode anexar e carregar código
compilado exatamente como se faz com código interpretado (no momento,
anexar e carregar código em Cython não é suportado no Notebook). A
compilação é realizada implicitamente sem que você tenha que fazer
qualquer coisa explicitamente. Veja
``$SAGE_ROOT/examples/programming/sagex/factorial.spyx`` para um
exemplo de uma implementação compilada para a função fatorial que usa
diretamente a biblioteca GMP em C. Experimente o seguinte, usando cd,
vá para o diretório ``$SAGE_ROOT/examples/programming/sagex/``, e
então faça o seguinte:

.. skip

::

    sage: load "factorial.spyx"
    ***************************************************
                    Recompiling factorial.spyx
    ***************************************************
    sage: factorial(50)
    30414093201713378043612608166064768844377641568960512000000000000
    sage: time n = factorial(10000)
    CPU times: user 0.03 s, sys: 0.00 s, total: 0.03 s
    Wall time: 0.03

Note que o Sage vai recompilar ``factorial.spyx`` se você encerrar e
reiniciar o Sage. A biblioteca compilada e compartilhada é armazenada
em ``$HOME/.sage/temp/hostname/pid/spyx``. Esses arquivos são
excluídos quando você encerra o Sage.

Nenhum pré-processamento é aplicado em arquivos spyx, por exemplo,
``1/3`` vai resultar em 0 em um arquivo spyx em vez do número racional
:math:`1/3`. Se ``foo`` é uma função da biblioteca Sage, para usá-la
em um arquivo spyx importe ``sage.all`` e use ``sage.all.foo``.

::

    import sage.all
    def foo(n):
        return sage.all.factorial(n)

Acessando Funções em C em Arquivos Separados
--------------------------------------------

É fácil também acessar funções em C definidas em arquivos \*.c
separados. Aqui vai um exemplo. Crie os arquivos ``test.c`` e
``test.spyx`` no mesmo diretório contendo:

Código C puro: ``test.c``

::

    int add_one(int n) {
      return n + 1;
    }

Código Cython: ``test.spyx``:

::

    cdef extern from "test.c":
        int add_one(int n)
    
    def test(n):
        return add_one(n)

Então o seguinte funciona:

.. skip

::

    sage: attach "test.spyx"
    Compiling (...)/test.spyx...
    sage: test(10)
    11

Se uma biblioteca ``foo`` adicional é necessária para compilar código
em C gerado a partir de um arquivo em Cython, adicione a linha ``clib
foo`` no arquivo fonte em Cython. De forma similar, um arquivo em C
adicional ``bar`` pode ser incluído na compilação declarando ``cfile
bar``.

.. _section-standalone:

Scripts Independentes em Python/Sage
====================================

O seguinte script em Sage fatora inteiros, polinômios, etc:

::

    #!/usr/bin/env sage
    
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: %s <n>" % sys.argv[0])
        print("Outputs the prime factorization of n.")
        sys.exit(1)
    
    print(factor(sage_eval(sys.argv[1])))

Para usar esse script, sua ``SAGE_ROOT`` precisa estar na sua variável
PATH. Se o script acima for chamado ``factor``, aqui está um exemplo
de como usá-lo:

::

    bash $ ./factor 2006
    2 * 17 * 59

Tipo de Dados
=============

Cada objeto em Sage possui um tipo bem definido. O Python possui
diversos tipos de dados, e a biblioteca do Sage adiciona ainda mais.
Os tipos de dados de Python incluem strings, listas, tuplas, inteiros
e floats, como ilustrado:

::

    sage: s = "sage"; type(s)
    <... 'str'>
    sage: s = 'sage'; type(s)      # you can use either single or double quotes
    <... 'str'>
    sage: s = [1,2,3,4]; type(s)
    <... 'list'>
    sage: s = (1,2,3,4); type(s)
    <... 'tuple'>
    sage: s = int(2006); type(s)
    <... 'int'>
    sage: s = float(2006); type(s)
    <... 'float'>

Além disso, o Sage acrescenta vários outros tipos. Por exemplo,
espaços vetoriais:

::

    sage: V = VectorSpace(QQ, 1000000); V
    Vector space of dimension 1000000 over Rational Field
    sage: type(V)
    <class 'sage.modules.free_module.FreeModule_ambient_field_with_category'>

Apenas certas funções podem ser aplicadas sobre ``V``. Em outros
softwares de matemática, essas seriam chamadas usando a notação
"funcional" ``foo(V,...)``. Em Sage, algumas funções estão anexadas ao
tipo (ou classe) de ``V``, e são chamadas usando uma sintaxe orientada
a objetos como em Java ou C++, por exemplo, ``V.foo()``. Isso ajuda a
manter o espaço de variáveis global sem milhares de funções, e permite
que várias funções diferentes com comportamento diferente possam ser
chamadas foo, sem a necessidade de usar um mecanismo de identificação
de tipos (ou casos) para decidir qual chamar. Além disso, se você
reutilizar o nome de uma função, essa função continua ainda disponível
(por exemplo, se você chamar algo ``zeta``, e então quiser calcular o
valor da função zeta de Riemann em 0.5, você continua podendo digitar
``s=.5; s.zeta()``).

::

    sage: zeta = -1
    sage: s=.5; s.zeta()     
    -1.46035450880959

Em alguns casos muito comuns, a notação funcional usual é também
suportada por conveniência e porque expressões matemáticas podem
parecer confusas usando a notação orientada a objetos. Aqui vão alguns
exemplos.

::

    sage: n = 2; n.sqrt()
    sqrt(2)
    sage: sqrt(2)
    sqrt(2)
    sage: V = VectorSpace(QQ,2)
    sage: V.basis()
        [(1, 0), (0, 1)]
    sage: basis(V)
        [(1, 0), (0, 1)]
    sage: M = MatrixSpace(GF(7), 2); M
    Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
    sage: A = M([1,2,3,4]); A
    [1 2]
    [3 4]
    sage: A.charpoly('x')
    x^2 + 2*x + 5
    sage: charpoly(A, 'x')
    x^2 + 2*x + 5

Para listar todas as funções para :math:`A`, use completamento tab.
Simplesmente digite ``A.``, então tecle ``[tab]`` no seu teclado, como
descrito em :ref:`section-tabcompletion`.

Listas, Tuplas e Sequências
===========================

O tipo de dados lista armazena elementos de um tipo arbitrário. Como
em C, C++, etc. (mas diferentemente da maioria dos sistemas de álgebra
computacional), os elementos da lista são indexados a partir do
:math:`0`:

::

    sage: v = [2, 3, 5, 'x', SymmetricGroup(3)]; v
    [2, 3, 5, 'x', Symmetric group of order 3! as a permutation group]
    sage: type(v)
    <... 'list'>
    sage: v[0]
    2
    sage: v[2]
    5

(Quando se indexa uma lista, é permitido que o índice não seja um int
do Python!) Um Inteiro do Sage (ou Racional, ou qualquer objeto que
possua um método ``__index__``) também ira funcionar.

::

    sage: v = [1,2,3]
    sage: v[2]
    3
    sage: n = 2      # Sage Integer
    sage: v[n]       # Perfectly OK!
    3
    sage: v[int(n)]  # Also OK.
    3

A função ``range`` cria uma lista de int's do Python (não Inteiros do
Sage):

::

    sage: list(range(1, 15))
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

Isso é útil quando se usa "list comprehensions" para construir listas:

::

    sage: L = [factor(n) for n in range(1, 15)]
    sage: L
    [1, 2, 3, 2^2, 5, 2 * 3, 7, 2^3, 3^2, 2 * 5, 11, 2^2 * 3, 13, 2 * 7]
    sage: L[12]
    13
    sage: type(L[12])
    <class 'sage.structure.factorization_integer.IntegerFactorization'>
    sage: [factor(n) for n in range(1, 15) if is_odd(n)]
    [1, 3, 5, 7, 3^2, 11, 13]

Para mais sobre como criar listas usando "list comprehensions", veja
[PyT]_.

Fatiamento de lista (list slicing) é um recurso fantástico. Se ``L`` é
uma lista, então ``L[m:n]`` retorna uma sub-lista de ``L`` obtida
começando do :math:`m`-ésimo elemento e terminando no
:math:`(n-1)`-ésimo elemento, como ilustrado abaixo.

::

    sage: L = [factor(n) for n in range(1, 20)]
    sage: L[4:9]
    [5, 2 * 3, 7, 2^3, 3^2]
    sage: L[:4]
    [1, 2, 3, 2^2]
    sage: L[14:4]
    []
    sage: L[14:]
    [3 * 5, 2^4, 17, 2 * 3^2, 19]

Tuplas são semelhantes à listas, exceto que elas são imutáveis: uma
vez criadas elas não podem ser alteradas.

::

    sage: v = (1,2,3,4); v
    (1, 2, 3, 4)
    sage: type(v)
    <... 'tuple'>
    sage: v[1] = 5
    Traceback (most recent call last):
    ...   
    TypeError: 'tuple' object does not support item assignment

Sequências são um terceiro tipo de dados do Sage semelhante a listas.
Diferentemente de listas e tuplas, Sequence não é um tipo de dados
nativo do Python. Por definição, uma sequência é mutável, mas usando o
método ``set_immutable`` da classe ``Sequence`` elas podem ser feitas
imutáveis, como mostra o exemplo a seguir. Todos os elementos da
sequência possuem um parente comum, chamado o universo da sequência.

::

    sage: v = Sequence([1,2,3,4/5])
    sage: v
    [1, 2, 3, 4/5]
    sage: type(v)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: type(v[1])
    <class 'sage.rings.rational.Rational'>
    sage: v.universe()
    Rational Field
    sage: v.is_immutable()
    False
    sage: v.set_immutable()
    sage: v[0] = 3
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Sequências são derivadas de listas e podem ser usadas em qualquer
lugar que listas são usadas.

::

    sage: v = Sequence([1,2,3,4/5])
    sage: isinstance(v, list)
    True
    sage: list(v)
    [1, 2, 3, 4/5]
    sage: type(list(v))
    <... 'list'>

Como um outro exemplo, bases para espaços vetoriais são sequências
imutáveis, pois é importante que elas não sejam modificadas.

::

    sage: V = QQ^3; B = V.basis(); B
    [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    sage: type(B)
    <class 'sage.structure.sequence.Sequence_generic'>
    sage: B[0] = B[1]
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.
    sage: B.universe()
    Vector space of dimension 3 over Rational Field

Dicionários
===========

Um dicionário (também chamado as vezes de lista associativa ou "hash
table") é um mapeamento de objetos em objetos arbitrários. (Exemplos
de objetos que admitem uma lista associativa são strings e números;
veja a documentação Python em https://docs.python.org/3/tutorial/index.html
para detalhes).

::

    sage: d = {1:5, 'sage':17, ZZ:GF(7)}
    sage: type(d)
    <... 'dict'>
    sage: list(d.keys())
    [1, 'sage', Integer Ring]
    sage: d['sage']
    17
    sage: d[ZZ]
    Finite Field of size 7
    sage: d[1]
    5

A terceira chave (key) ilustra como os índices de um dicionário podem
ser complicados, por exemplo, um anel de inteiros.

Você pode transformar o dicionário acima em uma lista com os mesmos
dados:

.. link

::

    sage: list(d.items())
    [(1, 5), ('sage', 17), (Integer Ring, Finite Field of size 7)]

É comum iterar sobre os pares em um dicionário:

:: 

    sage: d = {2:4, 3:9, 4:16}
    sage: [a*b for a, b in d.items()]
    [8, 27, 64]

Um dicionário não possui ordem, como o exemplo acima mostra.

Conjuntos
=========

O Python possui um tipo de conjuntos (set) nativo. O principal recurso
que ele oferece é a rápida verificação se um objeto está ou não em um
conjunto, juntamente com as operações comuns em conjuntos.

::

    sage: X = set([1,19,'a']);   Y = set([1,1,1, 2/3])
    sage: X   # random
    {1, 19, 'a'}
    sage: Y   # random
    {2/3, 1}
    sage: 'a' in X
    True
    sage: 'a' in Y
    False
    sage: X.intersection(Y)
    {1}

O Sage também possui o seu próprio tipo de dados para conjuntos que é
(em alguns casos) implementado usando o tipo nativo do Python, mas
possuir algumas funcionalidades adicionais. Crie um conjunto em Sage
usando ``Set(...)``. Por exemplo,

::

    sage: X = Set([1,19,'a']);   Y = Set([1,1,1, 2/3])
    sage: X   # random
    {'a', 1, 19}
    sage: Y   # random
    {1, 2/3}
    sage: X.intersection(Y)
    {1}
    sage: print(latex(Y))
    \left\{1, \frac{2}{3}\right\}
    sage: Set(ZZ)
    Set of elements of Integer Ring

Iteradores
==========

Iteradores foram adicionados recentemente ao Python e são
particularmente úteis em aplicações matemáticas. Aqui estão vários
exemplos; veja [PyT]_ para mais detalhes. Vamos criar um iterador
sobre o quadrados dos números inteiros até :math:`10000000`.

::

    sage: v = (n^2 for n in range(10000000))
    sage: next(v)
    0
    sage: next(v)
    1
    sage: next(v)
    4

Criamos agora um iterador sobre os primos da forma :math:`4p+1` com
:math:`p` também primo, e observamos os primeiros valores.

::

    sage: w = (4*p + 1 for p in Primes() if is_prime(4*p+1))
    sage: w         # in the next line, 0xb0853d6c is a random 0x number
    <generator object at 0xb0853d6c>
    sage: next(w)
    13
    sage: next(w)
    29
    sage: next(w)
    53

Certos anéis, por exemplo, corpos finitos e os inteiros possuem
iteradores associados a eles:

::

    sage: [x for x in GF(7)]
    [0, 1, 2, 3, 4, 5, 6]
    sage: W = ((x,y) for x in ZZ for y in ZZ)
    sage: next(W)
    (0, 0)
    sage: next(W)
    (0, 1)
    sage: next(W)
    (0, -1)

Laços, Funções, Enunciados de Controle e Comparações
====================================================

Nós já vimos alguns exemplos de alguns usos comuns de laços (loops)
``for``. Em Python, um laço ``for`` possui uma estrutura tabulada, tal
como

::

    >>> for i in range(5):
    ...      print(i)       
    ...
    0
    1
    2
    3
    4

Note que os dois pontos no final do enunciado (não existe "do" ou "od"
como no GAP ou Maple), e a identação antes dos comandos dentro do
laço, isto é, ``print(i)``. A tabulação é importante. No Sage, a
tabulação é automaticamente adicionada quando você digita ``enter``
após ":", como ilustrado abaixo.

::

    sage: for i in range(5):
    ....:     print(i)  # now hit enter twice
    0
    1
    2
    3
    4

O símbolo ``=`` é usado para atribuição.
O símbolo ``==`` é usado para verificar igualdade:

::

    sage: for i in range(15):
    ....:     if gcd(i,15) == 1:
    ....:         print(i)
    1
    2
    4
    7
    8
    11
    13
    14

Tenha em mente como a tabulação determina a estrutura de blocos para
enunciados ``if``, ``for``, e ``while``:

::

    sage: def legendre(a,p):
    ....:     is_sqr_modp=-1
    ....:     for i in range(p):
    ....:         if a % p == i^2 % p:
    ....:             is_sqr_modp=1
    ....:     return is_sqr_modp
             
    sage: legendre(2,7)
    1
    sage: legendre(3,7)
    -1

Obviamente essa não é uma implementação eficiente do símbolo de
Legendre! O objetivo é ilustrar vários aspectos da programação em
Python/Sage. A função {kronecker}, que já vem com o Sage, calcula o
símbolo de Legendre eficientemente usando uma biblioteca em C do PARI.

Finalmente, observamos que comparações, tais como ``==``, ``!=``,
``<=``, ``>=``, ``>``, ``<``, entre números irão automaticamente
converter ambos os números para o mesmo tipo, se possível:

::

    sage: 2 < 3.1; 3.1 <= 1
    True
    False
    sage: 2/3 < 3/2;   3/2 < 3/1
    True
    True

Use bool para desigualdades simbólicas:

::

    sage: x < x + 1
    x < x + 1
    sage: bool(x < x + 1)
    True

Quando se compara objetos de tipos diferentes no Sage, na maior parte
dos casos o Sage tenta encontrar uma coação canônica para ambos os
objetos em um parente comum (veja :ref:`section-coercion` para mais
detalhes). Se isso for bem sucedido, a comparação é realizada entre os
objetos que foram coagidos; se não for bem sucedido, os objetos são
considerados diferentes. Para testar se duas variáveis fazem
referência ao mesmo objeto use ``is``. Como se vê no próximo exemplo,
o int ``1`` do Python é único, mas o Inteiro ``1`` do Sage não é.

::

    sage: 1 is 2/2
    False
    sage: 1 is 1
    False
    sage: 1 == 2/2
    True

Nas duas linhas seguintes, a primeira igualdade é falsa (``False``)
porque não existe um morfismo canônico :math:`QQ\to \GF{5}`, logo
não há uma forma de comparar o :math:`1` em :math:`\GF{5}` com o
:math:`1 \in \QQ`. Em contraste, existe um mapa canônico entre
:math:`\ZZ \to \GF{5}`, logo a segunda comparação é verdadeira
(``True``)

::

    sage: GF(5)(1) == QQ(1); QQ(1) == GF(5)(1)
    False
    False
    sage: GF(5)(1) == ZZ(1); ZZ(1) == GF(5)(1)
    True
    True
    sage: ZZ(1) == QQ(1)
    True

ATENÇÃO: Comparação no Sage é mais restritiva do que no Magma, o qual
declara :math:`1 \in \GF{5}` igual a :math:`1 \in \QQ`.

::

    sage: magma('GF(5)!1 eq Rationals()!1')            # optional magma required
    true

Otimização (Profiling)
======================

Autor desta seção: Martin Albrecht (https://martinralbrecht.wordpress.com/)

    "Premature optimization is the root of all evil." - Donald Knuth

As vezes é útil procurar por gargalos em programas para entender quais
partes gastam maior tempo computacional; isso pode dar uma boa ideia
sobre quais partes otimizar. Python e portanto Sage fornecem várias
opções de "profiling" (esse é o nome que se dá ao processo de
otimização).

O mais simples de usar é o comando ``prun`` na linha de comando
interativa. Ele retorna um sumário sobre o tempo computacional
utilizado por cada função. Para analisar (a atualmente lenta! -- na
versão 1.0) multiplicação de matrizes sobre corpos finitos, por
exemplo, faça o seguinte:

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])

.. skip

::

    sage: %prun B = A*A
           32893 function calls in 1.100 CPU seconds
    
    Ordered by: internal time
    
    ncalls tottime percall cumtime percall filename:lineno(function)
     12127  0.160   0.000   0.160  0.000 :0(isinstance)
      2000  0.150   0.000   0.280  0.000 matrix.py:2235(__getitem__)
      1000  0.120   0.000   0.370  0.000 finite_field_element.py:392(__mul__)
      1903  0.120   0.000   0.200  0.000 finite_field_element.py:47(__init__)
      1900  0.090   0.000   0.220  0.000 finite_field_element.py:376(__compat)
       900  0.080   0.000   0.260  0.000 finite_field_element.py:380(__add__)
         1  0.070   0.070   1.100  1.100 matrix.py:864(__mul__)
      2105  0.070   0.000   0.070  0.000 matrix.py:282(ncols)
      ...

Aqui ``ncalls`` é o números de chamadas, ``tottime`` é o tempo total
gasto por uma determinada função (excluíndo o tempo gasto em chamadas
de subfunções), ``percall`` é o quociente de ``tottime`` dividido por
``ncalls``. ``cumtime`` é o tempo total gasto nessa e em todas as
subfunções (isto é, desde o início até o término da execução da
função), ``percall`` é o quociente de ``cumtime`` dividido pelas
chamadas primitivas, e ``filename:lineno(function)`` fornece os dados
respectivos para cada função. A regra prática aqui é: Quanto mais no
topo uma função aparece nessa lista, mais custo computacional ela
acarreta. Logo é mais interessante para ser optimizada.

Como usual, ``prun?`` fornece detalhes sobre como usar o "profiler" e
como entender a saída de dados.

A saída de dados pode ser escrita em um objeto para permitir uma
análise mais detalhada:

.. skip

::

    sage: %prun -r A*A
    sage: stats = _
    sage: stats?

Note: digitando ``stats = prun -r A\*A`` obtém-se um erro de sintaxe
porque prun é um comando do IPython, não uma função comum.

Para uma representação gráfica dos dados do "profiling", você pode
usar o "hotspot profiler", um pequeno script chamado
``hotshot2cachetree`` e o programa ``kcachegrind`` (apenas no Unix). O
mesmo exemplo agora com o "hotspot profiler":

.. skip

::

    sage: k,a = GF(2**8, 'a').objgen()
    sage: A = Matrix(k,10,10,[k.random_element() for _ in range(10*10)])
    sage: import hotshot
    sage: filename = "pythongrind.prof"
    sage: prof = hotshot.Profile(filename, lineevents=1)

.. skip

::

    sage: prof.run("A*A")
    <hotshot.Profile instance at 0x414c11ec>
    sage: prof.close()

Isso resulta em um arquivo ``pythongrind.prof`` no diretório de
trabalho atual. Ele pode ser convertido para o formato cachegrind para
visualização.

Em uma linha de comando do sistema, digite

.. skip

::

    hotshot2calltree -o cachegrind.out.42 pythongrind.prof

O arquivo de saída ``cachegrind.out.42`` pode ser examinado com
``kcachegrind``. Note que a convenção de nomes ``cachegrind.out.XX``
precisa ser obedecida.
