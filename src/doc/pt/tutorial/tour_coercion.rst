.. -*- coding: utf-8 -*-

.. _section-coercion:

============================
Famílias, Conversão e Coação
============================

Esta seção pode parecer mais técnica do que as anteriores, mas
acreditamos que é importante entender o significado de famílias e
coação de modo a usar anéis e outras estruturas algébricas no Sage de
forma efetiva e eficiente.

Note que vamos explicar algumas noções, mas não vamos mostrar aqui
como implementá-las. Um tutorial voltado à implementação está
disponível (em inglês) como um 
`tutorial temática <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_.

Elementos
---------

Caso se queira implementar um anel em Python, uma primeira aproximação
seria criar uma classe para os elementos ``X`` do anel e adicionar os
requeridos métodos (com underscores duplos) ``__add__``, ``__sub``,
``__mul__``, obviamente garantindo que os axiomas de anel são
verificados.

Como o Python é uma linguagem de tipagem forte (ainda que de tipagem
dinâmica), poderia-se, pelo menos a princípio, esperar-se que fosse
implementado em Python uma classe para cada anel. No final das contas,
o Python contém um tipo ``<int>`` para os inteiros, um tipo
``<float>`` para os reais, e assim por diante. Mas essa estratégia
logo encontra uma limitação: Existe um número infinito de anéis, e não
se pode implementar um número infinito de classes.

Em vez disso, poderia-se criar uma hierarquia de classes projetada
para implementar elementos de estruturas algébricas ubíquas, tais como
grupos, anéis, anéis comutativos, corpos, álgebras, e assim por
diante.

Mas isso significa que elementos de anéis bastante diferentes podem
ter o mesmo tipo.

::

    sage: P.<x,y> = GF(3)[]
    sage: Q.<a,b> = GF(4,'z')[]
    sage: type(x)==type(a)
    True

Por outro lado, poderia-se ter também classes diferentes em Python
fornecendo implementações diferentes da mesma estrutura matemática
(por exemplo, matrizes densas versus matrizes esparsas).

::

    sage: P.<a> = PolynomialRing(ZZ)
    sage: Q.<b> = PolynomialRing(ZZ, sparse=True)
    sage: R.<c> = PolynomialRing(ZZ, implementation='NTL')
    sage: type(a) == type(b)
    False
    sage: type(a) == type(c)
    False
    sage: type(b) == type(c)
    False

Isso apresenta dois problemas: Por um lado, se tivéssemos elementos
que são duas instancias da mesma classe, então poderia-se esperar que
o método ``__add__`` dessas classes permitisse somá-los; mas não
se deseja isso, se os elementos pertencem a anéis bastante diferentes.
Por outro lado, se possui-se elementos que pertencem a implementações
diferentes do mesmo anel, então gostaria-se de somá-los, mas isso não
pode ser feito diretamente se eles pertencem a classes diferentes em
Python.

A solução para esses problemas é chamada coação e será explicada a
seguir.

Todavia, é essencial que cada elemento saiba a qual pertence. Isso
está disponível através método ``parent()``:

.. link

::

    sage: a.parent(); b.parent(); c.parent()
    Univariate Polynomial Ring in a over Integer Ring
    Sparse Univariate Polynomial Ring in b over Integer Ring
    Univariate Polynomial Ring in c over Integer Ring (using NTL)


Famílias e Categorias
---------------------

De forma similar à hierarquia de classes em Python voltada para
elementos de estruturas algébricas, o Sage também fornece classes para
as estruturas algébricas que contém esses elementos. Estruturas
contendo elementos são chamadas "estruturas parente" no Sage, e existe
uma classe básica para elas. Paralelamente à hierarquia de noções
matemáticas, tem-se uma hierarquia de classes, a saber, para
conjuntos, anéis, corpos e assim por diante:

::

    sage: isinstance(QQ,Field)
    True
    sage: isinstance(QQ, Ring)
    True
    sage: isinstance(ZZ,Field)
    False
    sage: isinstance(ZZ, Ring)
    True

Em álgebra, objetos que compartilham o mesmo tipo de estruturas
algébricas são agrupados nas assim chamadas "categorias". Logo, existe
uma analogia aproximada entre a hierarquia de classes em Sage e a
hierarquia de categorias. Todavia, essa analogia de classes em Python
e categorias não deve ser enfatizada demais. No final das contas,
categorias matemáticas também são implementadas no Sage:

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

Enquanto a hierarquia de classes no Sage é centrada nos detalhes de
implementação, a construção de categorias em Sage é mais centrada
na estrutura matemática. É possível implementar métodos e testes
gerais independentemente de uma implementação específica nas
categorias.

Estruturas da mesma família em Sage são supostamente objetos únicos em
Python. Por exemplo, uma vez que um anel de polinômios sobre um certo anel
base e com uma certa lista de geradores é criada, o resultado é arquivado:

::

    sage: RR['x','y'] is RR['x','y']
    True


Tipos versus Parentes
---------------------

O tipo ``RingElement`` não deve ser confundido com a noção matemática
de elemento de anel; por razões práticas, as vezes um objeto é uma
instancia de ``RingElement`` embora ele não pertence a um anel:

::

    sage: cristovao = ZZ(1492)
    sage: isinstance(cristovao, RingElement)
    True

Enquanto *famílias* são únicas, elementos iguais de uma família em Sage
não são necessariamente idênticos. Isso contrasta com o comportamento
do Python para alguns (embora não todos) inteiros:

::

    sage: int(1) is int(1) # Python int
    True
    sage: int(-15) is int(-15)
    False
    sage: 1 is 1           # Sage Integer
    False

É importante observar que elementos de anéis diferentes em geral não
podem ser distinguidos pelos seus tipos, mas sim por sua família:

::

    sage: a = GF(2)(1)
    sage: b = GF(5)(1)
    sage: type(a) is type(b)
    True
    sage: parent(a)
    Finite Field of size 2
    sage: parent(b)
    Finite Field of size 5

Logo, de um ponto de vista algébrico, **o parente de um elemento é
mais importante do que seu tipo.**

Conversão versus Coação
-----------------------

Em alguns casos é possível converter um elemento de uma estrutura
parente em um elemento de uma outra estrutura parente. Tal conversão
pode ser tanto explícita como implícita (essa é chamada *coação*).

O leitor pode conhecer as noções de *conversão de tipo* e *coação de
tipo* como na linguagem C, por exemplo. Existem noções de *conversão*
e *coação* em Sage também. Mas as noções em Sage são centradas em
*família*, não em tipos. Então, por favor não confunda conversão de
tipo em C com conversão em Sage!

Aqui se encontra uma breve apresentação. Para uma descrição detalhada
e informações sobre a implementação, referimos à seção sobre coação no
manual de referência e para o `tutorial
<http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_.

Existem duas possibilidades extremas com respeito à possibilidade de
fazer aritmética com elementos de *anéis diferentes*:

* Anéis diferentes são mundos diferentes, e não faz nenhum sentido
  somar ou multiplicar elementos de anéis diferentes; mesmo ``1 +
  1/2`` não faz sentido, pois o primeiro somando é um inteiro e o
  segundo um racional.

Ou

* Se um elemento ``r1`` de uma aner ``R1`` pode de alguma forma ser
  interpretado em um outro anel ``R2``, então todas as operações
  aritméticas envolvendo ``r1`` e qualquer elemento de ``R2`` são
  permitidas. O elemento neutro da multiplicação existe em todos os
  corpos e em vários anéis, e eles devem ser todos iguais.

O Sage faz uma concessão. Se ``P1`` e ``P2`` são estruturas da mesma família
e ``p1`` é um elemento de ``P1``, então o usuário pode explicitamente
perguntar por uma interpretação de ``p1`` em ``P2``. Isso pode não fazer
sentido em todos os casos ou não estar definido para todos os elementos de
``P1``, e fica a cargo do usuário assegurar que isso faz sentido. Nos
referimos a isso como **conversão**:

::

    sage: a = GF(2)(1)
    sage: b = GF(5)(1)
    sage: GF(5)(a) == b
    True
    sage: GF(2)(b) == a
    True

Todavia, uma conversão *implícita* (ou automática) ocorrerá apenas se
puder ser feita *completamente* e *consistentemente*. Rigor matemático
é essencial nesse ponto.

Uma tal conversão implícita é chamada **coação**. Se coação for
definida, então deve coincidir com conversão. Duas condições devem ser
satisfeitas para uma coação ser definida:

#. Uma coação de ``P1`` para ``P2`` deve ser dada por uma estrutura
   que preserva mapeamentos (por exemplo, um homomorfismo de anéis).
   Não é suficiente que *alguns* elementos de ``P1`` possam ser
   mapeados em ``P2``, e o mapa deve respeitar a estrutura algébrica
   de ``P1``.
#. A escolha desses mapas de coação deve ser consistente: Se ``P3`` é
   uma terceira estrutura parente, então a composição da coação
   adotada de ``P1`` para ``P2`` com a coação de ``P2`` para ``P3``
   deve coincidir com a coação adotada de ``P1`` para ``P3``. Em
   particular, se existir uma coação de ``P1`` para ``P2`` e ``P2``
   para ``P1``, a composição deve ser o mapa identidade em ``P1``.

Logo, embora é possível converter cada elemento de ``GF(2)`` para
``GF(5)``, não há coação, pois não existe homomorfismo de anel entre
``GF(2)`` e ``GF(5)``.

O segundo aspecto - consistência - é um pouco mais difícil de
explicar. Vamos ilustrá-lo usando anéis de polinômios em mais de uma
variável. Em aplicações, certamente faz mais sentido ter coações que
preservam nomes. Então temos:

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

Se não existir homomorfismo de anel que preserve nomes, coação não é
definida. Todavia, conversão pode ainda ser possível, a saber,
mapeando geradores de anel de acordo com sua posição da lista de
geradores:

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

Mas essas conversões que preservam a posição não se qualificam como
coação: Compondo um mapa que preserva nomes de ``ZZ['x','y']`` para
``ZZ['y','x']``, com um mapa que preserva nomes de ``ZZ['y','x']``
para ``ZZ['a','b']``, resultaria em um mapa que não preserva nomes nem
posição, violando a consistência.

Se houver coação, ela será usada para comparar elementos de anéis
diferentes ou fazer aritmética. Isso é frequentemente conveniente, mas
o usuário deve estar ciente que estender a relação ``==`` além das
fronteiras de famílias diferentes pode facilmente resultar em 
problemas. Por exemplo, enquanto ``==`` é supostamente uma relação de
equivalência sobre os elementos de *um* anel, isso não é
necessariamente o caso se anéis *diferentes* estão envolvidos. Por
exemplo, ``1`` em ``ZZ`` e em um corpo finito são considerados iguais,
pois existe uma coação canônica dos inteiros em qualquer corpo finito.
Todavia, em geral não existe coação entre dois corpos finitos
diferentes. Portanto temos

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

Similarmente,

.. link

::

    sage: R3(R1.1) == R3.1
    True
    sage: R1.1 == R3.1
    False
    sage: R1.1 != R3.1
    True

Uma outra consequência da condição de consistência é que coação pode
apenas ir de anéis exatos (por exemplo, os racionais ``QQ``) para
anéis não-exatos (por exemplo, os números reais com uma precisão fixa
``RR``), mas não na outra direção. A razão é que a composição da
coação de ``QQ`` em ``RR`` com a conversão de ``RR`` para ``QQ``
deveria ser a identidade em ``QQ``. Mas isso é impossível, pois alguns
números racionais distintos podem ser tratados como iguais em ``RR``,
como no seguinte exemplo:

::

    sage: RR(1/10^200+1/10^100) == RR(1/10^100)
    True
    sage: 1/10^200+1/10^100 == 1/10^100
    False

Quando se compara elementos de duas famílias ``P1`` e ``P2``, é
possível que não haja coação entre os dois anéis, mas existe uma
escolha canônica de um parente ``P3`` de modo que tanto ``P1`` como
``P2`` são coagidos em ``P3``. Nesse caso, coação vai ocorrer também.
Um caso de uso típico é na soma de um número racional com um polinômio
com coeficientes inteiros, resultando em um polinômio com coeficientes
racionais:

::

    sage: P1.<x> = ZZ[]
    sage: p = 2*x+3
    sage: q = 1/2
    sage: parent(p)
    Univariate Polynomial Ring in x over Integer Ring
    sage: parent(p+q)
    Univariate Polynomial Ring in x over Rational Field

Note que a princípio o resultado deveria também fazer sentido no
corpo de frações de ``ZZ['x']``. Todavia, o Sage tenta escolher um
parente *canônico* comum que parece ser o mais natural (``QQ['x']`` no
nosso exemplo). Se várias famílias potencialmente comuns parecem
igualmente naturais, o Sage *não* vai escolher um deles
aleatoriamente. Os mecanismos sobre os quais essa escolha se baseia é
explicado em um
`tutorial <http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html>`_

Nenhuma coação para um parente comum vai ocorrer no seguinte exemplo:

::

    sage: R.<x> = QQ[]
    sage: S.<y> = QQ[]
    sage: x+y
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +: 'Univariate Polynomial Ring in x over Rational Field' and 'Univariate Polynomial Ring in y over Rational Field'

A razão é que o Sage não escolhe um dos potenciais candidatos
``QQ['x']['y']``, ``QQ['y']['x']``, ``QQ['x','y']`` ou
``QQ['y','x']``, porque todas essas estruturas combinadas em pares
diferentes parecem ser de famílias comuns naturais, e não existe escolha
canônica aparente.
