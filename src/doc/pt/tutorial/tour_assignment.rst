Atribuição, Igualdade, e Aritmética
===================================

Com pequenas exceções, o Sage utiliza a linguagem de programação
Python, logo a maioria dos livros de introdução ao Python vão ajudá-lo
a aprender Sage.

O Sage usa ``=`` para atribuição, e usa ``==``, ``<=``, ``>=``, ``<``
e ``>`` para comparação:

::

    sage: a = 5
    sage: a
    5
    sage: 2 == 2
    True
    sage: 2 == 3
    False
    sage: 2 < 3
    True
    sage: a == 5
    True

O Sage fornece todas as operações matemáticas básicas:

::

    sage: 2**3    #  ** means exponent
    8
    sage: 2^3     #  ^ is a synonym for ** (unlike in Python)
    8
    sage: 10 % 3  #  for integer arguments, % means mod, i.e., remainder
    1
    sage: 10/4
    5/2
    sage: 10//4   #  for integer arguments, // returns the integer quotient
    2
    sage: 4 * (10 // 4) + 10 % 4 == 10
    True
    sage: 3^2*4 + 2%5
    38

O cálculo de uma expressão como ``3^2*4 + 2%5`` depende da ordem em
que as operações são aplicadas; isso é especificado na "tabela de
precedência" em :ref:`section-precedence`.

O Sage também fornece várias funções matemáticas básicas; aqui estão
apenas alguns exemplos:

::

    sage: sqrt(3.4)
    1.84390889145858 
    sage: sin(5.135)
    -0.912021158525540 
    sage: sin(pi/3)
    1/2*sqrt(3)

Como o último exemplo mostra, algumas expressões matemáticas retornam
valores 'exatos' em vez de aproximações numéricas. Para obter uma
aproximação numérica, use a função ``n`` ou o método ``n`` (ambos
possuem um nome longo, ``numerical_approx``, e a função ``N`` é o
mesma que ``n``). Essas funções aceitam o argumento opcional
``prec``, que é o número de bits de precisão requisitado, e
``digits``, que é o número de dígitos decimais de precisão
requisitado; o padrão é 53 bits de precisão.

::

    sage: exp(2)
    e^2
    sage: n(exp(2))
    7.38905609893065
    sage: sqrt(pi).numerical_approx()
    1.77245385090552
    sage: sin(10).n(digits=5)
    -0.54402
    sage: N(sin(10),digits=10)
    -0.5440211109 
    sage: numerical_approx(pi, prec=200)
    3.1415926535897932384626433832795028841971693993751058209749

O Python é uma linguagem de tipagem dinâmica, portanto o valor
referido por cada variável possui um tipo associado a ele, mas uma
variável pode possuir valores de qualquer tipo em determinado escopo:

::

    sage: a = 5   # a is an integer
    sage: type(a)
    <class 'sage.rings.integer.Integer'>
    sage: a = 5/3  # now a is a rational number
    sage: type(a)
    <class 'sage.rings.rational.Rational'>
    sage: a = 'hello'  # now a is a string
    sage: type(a)
    <... 'str'>

A linguagem de programação C, que é de tipagem estática , é muito
diferente; uma variável que foi declarada como int pode apenas
armazenar um int em seu escopo.
