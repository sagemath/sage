.. _chapter-help:

Obtendo ajuda
=============

O Sage possui vasta documentação, acessível digitando o nome de uma
função ou constante (por exemplo), seguido pelo ponto de interrogação:

.. skip

::

    sage: tan?
    Type:        <class 'sage.calculus.calculus.Function_tan'>
    Definition:  tan( [noargspec] )
    Docstring: 
    
        The tangent function
    
        EXAMPLES:
            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
    sage: log2?
    Type:        <class 'sage.functions.constants.Log2'>
    Definition:  log2( [noargspec] )
    Docstring: 
    
        The natural logarithm of the real number 2.
        
        EXAMPLES:
            sage: log2
            log2
            sage: float(log2)
            0.69314718055994529
            sage: RR(log2)
            0.693147180559945
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(log2)
            0.69314718055994530941723212145817656807550013436025525412068
            sage: l = (1-log2)/(1+log2); l
            (1 - log(2))/(log(2) + 1)
            sage: R(l)
            0.18123221829928249948761381864650311423330609774776013488056
            sage: maxima(log2)
            log(2)
            sage: maxima(log2).float()
            .6931471805599453
            sage: gp(log2)
            0.6931471805599453094172321215             # 32-bit
            0.69314718055994530941723212145817656807   # 64-bit
    sage: sudoku?
    File:        sage/local/lib/python2.5/site-packages/sage/games/sudoku.py
    Type:        <... 'function'>
    Definition:  sudoku(A)
    Docstring: 
    
        Solve the 9x9 Sudoku puzzle defined by the matrix A.
    
        EXAMPLE:
            sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0,
        0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0, 0,0,0, 2,0,8, 0,0,0,
        0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2,
        0,0,0, 4,9,0, 0,5,0, 0,0,3])
            sage: A
            [5 0 0 0 8 0 0 4 9]
            [0 0 0 5 0 0 0 3 0]
            [0 6 7 3 0 0 0 0 1]
            [1 5 0 0 0 0 0 0 0]
            [0 0 0 2 0 8 0 0 0]
            [0 0 0 0 0 0 0 1 8]
            [7 0 0 0 0 4 1 5 0]
            [0 3 0 0 0 2 0 0 0]
            [4 9 0 0 5 0 0 0 3]
            sage: sudoku(A)
            [5 1 3 6 8 7 2 4 9]
            [8 4 9 5 2 1 6 3 7]
            [2 6 7 3 4 9 5 8 1]
            [1 5 8 4 6 3 9 7 2]
            [9 7 4 2 1 8 3 6 5]
            [3 2 6 7 9 5 4 1 8]
            [7 8 2 9 3 4 1 5 6]
            [6 3 5 1 7 2 8 9 4]
            [4 9 1 8 5 6 7 2 3]

O Sage também fornece completamento tab: digite as primeiras letras de
uma função e então pressione a tecla tab. Por exemplo, se você digitar
``ta`` seguido de ``TAB``, o Sage vai imprimir ``tachyon, tan, tanh,
taylor``. Essa é uma boa forma de encontrar nomes de funções e outras
estruturas no Sage.


.. _section-functions:

Funções, Tabulação, e Contagem
===============================

Para definir uma nova função no Sage, use o comando ``def`` e dois
pontos após a lista de nomes das variáveis. Por exemplo:

::

    sage: def is_even(n):
    ....:     return n % 2 == 0
    ....:
    sage: is_even(2)
    True
    sage: is_even(3)
    False

Observação: Dependendo da versão do tutorial que você está lendo,
você pode ver três pontos ``....:`` na segunda linha desse exemplo. Não
digite esses pontos; eles são apenas para enfatizar que o código está
tabulado. Se for esse o caso, pressione [Enter] uma vez após o fim do
bloco de código para inserir uma linha em branco e concluir a
definição da função.

Você não especifica o tipo de dado de nenhum dos argumentos da função.
É possível especificar argumentos múltiplos, cada um dos quais pode
ter um valor opcional padrão. Por exemplo, a função abaixo usa o valor
padrão ``divisor=2`` se ``divisor`` não é especificado.

::

    sage: def is_divisible_by(number, divisor=2):
    ....:     return number%divisor == 0
    sage: is_divisible_by(6,2)
    True
    sage: is_divisible_by(6)
    True
    sage: is_divisible_by(6, 5)
    False

Você também pode especificar explicitamente um ou mais argumentos
quando evocar uma função; se você especificar os argumentos
explicitamente, você pode fazê-lo em qualquer ordem:

.. link

::

    sage: is_divisible_by(6, divisor=5)
    False
    sage: is_divisible_by(divisor=2, number=6)
    True

Em Python, blocos de código não são indicados por colchetes ou blocos
de início e fim, como em outras linguagens. Em vez disso, blocos de
código são indicados por tabulação, que devem estar alinhadas
exatamente. Por exemplo, o seguinte código possui um erro de sintaxe
porque o comando ``return`` não possui a mesma tabulação da linha que
inicia o seu bloco de código.

.. skip

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:    return v
    Syntax Error:
           return v

Se você corrigir a tabulação, a função fica correta:

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:     return v
    sage: even(10)
    [4, 6, 8]

Não é necessário inserir ponto-e-vírgula no final da linha. Todavia,
você pode inserir múltiplos comandos em uma mesma linha separados por
ponto-e-vírgula:

::

    sage: a = 5; b = a + 3; c = b^2; c
    64

Se você quiser que uma única linha de comando seja escrita em mais de
uma linha, use ``\`` para quebrar a linha:

::

    sage: 2 + \
    ....:    3
    5

Em Sage, a contagem é feita iterando sobre um intervalo de inteiros.
Por exemplo, a primeira linha abaixo é equivalente a ``for(i=0; i<3;
i++)`` em C++ ou Java:

::

    sage: for i in range(3):
    ....:     print(i)
    0
    1
    2

A primeira linha abaixo é equivalente a ``for(i=2; i<5; i++)``.

::

    sage: for i in range(2,5):
    ....:     print(i)
    2
    3
    4

O Terceiro argumento controla o passo. O comando abaixo é equivalente
a ``for(i=1; i<6; i+=2)``.

::

    sage: for i in range(1,6,2):
    ....:     print(i)
    1
    3
    5

Frequentemente deseja-se criar uma tabela para visualizar resultados
calculados com o Sage. Uma forma fácil de fazer isso é utilizando
formatação de strings. Abaixo, criamos três colunas cada uma com
largura exatamente 6, e fazemos uma tabela com quadrados e cubos de
alguns números.

::

    sage: for i in range(5):
    ....:     print('%6s %6s %6s' % (i, i^2, i^3))
         0      0      0
         1      1      1
         2      4      8
         3      9     27
         4     16     64

A estrutura de dados mais básica em Sage é a lista, que é -- como o
nome sugere -- simplesmente uma lista de objetos arbitrários. Por
exemplo, o comando ``range`` que usamos acima cria uma lista:

::

    sage: list(range(2,10))
    [2, 3, 4, 5, 6, 7, 8, 9]

Abaixo segue uma lista mais complicada:

::

    sage: v = [1, "hello", 2/3, sin(x^3)]
    sage: v
    [1, 'hello', 2/3, sin(x^3)]

Listas são indexadas começando do 0, como em várias linguagens de
programação.

.. link

::

    sage: v[0]
    1
    sage: v[3]
    sin(x^3)

Use ``len(v)`` para obter o comprimento de ``v``, use
``v.append(obj)`` para inserir um novo objeto no final de ``v``, e use
``del v[i]`` para remover o :math:`i`-ésimo elemento de ``v``:

.. link

::

    sage: len(v)
    4
    sage: v.append(1.5)
    sage: v
    [1, 'hello', 2/3, sin(x^3), 1.50000000000000]
    sage: del v[1]
    sage: v
    [1, 2/3, sin(x^3), 1.50000000000000]

Outra importante estrutura de dados é o dicionário (ou lista
associativa). Ele funciona como uma lista, exceto que pode ser
indexado por vários tipos de objeto (os índices devem ser imutáveis):

::

    sage: d = {'hi':-2,  3/8:pi,   e:pi}
    sage: d['hi']
    -2
    sage: d[e]
    pi

Você pode também definir novos tipos de dados usando classes.
Encapsular objetos matemáticos usando classes é uma técnica poderosa
que pode ajudar a simplificar e organizar os seus programas em Sage.
Abaixo, definimos uma nova classe que representa a lista de inteiros
pares positivos até *n*; essa classe é derivada do tipo ``list``.

::

    sage: class Evens(list):
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:         list.__init__(self, range(2, n+1, 2))
    ....:     def __repr__(self):
    ....:         return "Even positive numbers up to n."

O método ``__init__`` é evocado para inicializar o objeto quando ele é
criado; o método ``__repr__`` imprime o objeto. Nós evocamos o
construtor ``__init__`` do tipo ``list`` na segunda linha do método
``__init__``. Criamos um objeto da classe ``Evens`` da seguinte forma:

.. link

::

    sage: e = Evens(10)
    sage: e
    Even positive numbers up to n.

Note que ``e`` imprime usando o método ``__repr__`` que nós
definimos. Para ver a lista de números, use a função ``list``:

.. link

::

    sage: list(e)
    [2, 4, 6, 8, 10]

Podemos também acessar o atributo ``n`` ou tratar ``e`` como uma
lista.

.. link

::

    sage: e.n
    10
    sage: e[2]
    6
