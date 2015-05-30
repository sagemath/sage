.. _section-poly:

Polinômios
==========

Nesta seção vamos ilustrar como criar e usar polinômios no Sage.


.. _section-univariate:

Polinômios em Uma Variável
--------------------------

Existem três formas de criar anéis de polinômios.

::

    sage: R = PolynomialRing(QQ, 't')
    sage: R
    Univariate Polynomial Ring in t over Rational Field

Esse comando cria um anel de polinômios e diz para o Sage usar a letra
't' para representar a variável indeterminada quando imprimir na tela.
Todavia, isso não define o símbolo ``t`` para uso no Sage, logo você
não pode usá-lo para definir um polinômio (como :math:`t^2+1`)
pertencente a ``R``.

Uma forma alternativa é

.. link

::

    sage: S = QQ['t']
    sage: S == R
    True

As mesmas observações com respeito a ``t`` valem também nesse caso.

Uma terceira e conveniente forma de definir polinômios é

::

    sage: R.<t> = PolynomialRing(QQ)

ou

::

    sage: R.<t> = QQ['t']

ou ainda

::

    sage: R.<t> = QQ[]

Isso tem o efeito colateral de definir a variável ``t`` como a
variável indeterminada do anel de polinômios, logo você pode
facilmente construir elementos de ``R`` da seguinte forma. (Note que
essa terceira alternativa é muito semelhante à notação usada em Magma,
e da mesma forma que no Magma ela pode ser usada para diversos tipos
de objetos.)

.. link

::

    sage: poly = (t+1) * (t+2); poly
    t^2 + 3*t + 2
    sage: poly in R
    True

Qualquer que seja o método usado para definir um anel de polinômios,
você pode recuperar a variável indeterminada como o :math:`0`-ésimo
gerador:

::

    sage: R = PolynomialRing(QQ, 't')
    sage: t = R.0
    sage: t in R
    True

Note que uma construção similar funciona com os números complexos: os
números complexos podem ser vistos como sendo gerados pelo símbolo
``i`` sobre os números reais; logo temos o seguinte:

::

    sage: CC
    Complex Field with 53 bits of precision
    sage: CC.0  # 0th generator of CC
    1.00000000000000*I

Para anel de polinômios, você pode obter tanto o anel como o seu
gerador, ou somente o gerador, no momento em que o anel for criado, da
seguinte forma:

::

    sage: R, t = QQ['t'].objgen()
    sage: t    = QQ['t'].gen()
    sage: R, t = objgen(QQ['t'])
    sage: t    = gen(QQ['t'])

Finalmente apresentamos um pouco de aritmética em :math:`\QQ[t]`.

::

    sage: R, t = QQ['t'].objgen()
    sage: f = 2*t^7 + 3*t^2 - 15/19
    sage: f^2
    4*t^14 + 12*t^9 - 60/19*t^7 + 9*t^4 - 90/19*t^2 + 225/361
    sage: cyclo = R.cyclotomic_polynomial(7); cyclo
    t^6 + t^5 + t^4 + t^3 + t^2 + t + 1
    sage: g = 7 * cyclo * t^5 * (t^5 + 10*t + 2)
    sage: g
    7*t^16 + 7*t^15 + 7*t^14 + 7*t^13 + 77*t^12 + 91*t^11 + 91*t^10 + 84*t^9 
           + 84*t^8 + 84*t^7 + 84*t^6 + 14*t^5
    sage: F = factor(g); F
    (7) * t^5 * (t^5 + 10*t + 2) * (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1)
    sage: F.unit()
    7
    sage: list(F)
    [(t, 5), (t^5 + 10*t + 2, 1), (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1, 1)]

Note que a fatorização corretamente leva em conta e armazena a parte
unitária.

Se você fosse usar, por exemplo, a função ``R.cyclotomic_polynomial``
intensamente para algum projeto de pesquisa, além de citar o Sage,
você deveria tentar descobrir qual componente do Sage é de fato usado
para calcular esses polinômios, e citá-lo também. Nesse caso, se você
digitar ``R.cyclotomic_polynomial??`` para ver o código fonte, você
irá facilmente ver uma linha ``f = pari.polcyclo(n)`` o que significa
que o PARI é usado para o cálculo dos polinômios ciclotrômicos. Cite o
PARI também no seu trabalho.

Dividindo dois polinômios cria-se um elemento do corpo de frações (o
qual o Sage cria automaticamente).

::

    sage: x = QQ['x'].0
    sage: f = x^3 + 1; g = x^2 - 17
    sage: h = f/g;  h
    (x^3 + 1)/(x^2 - 17)
    sage: h.parent()
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

Usando-se a série de Laurent, pode-se calcular a expansão em série no
corpo de frações de ``QQ[x]``:

::

    sage: R.<x> = LaurentSeriesRing(QQ); R
    Laurent Series Ring in x over Rational Field
    sage: 1/(1-x) + O(x^10)
    1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)

Se nomearmos a variável de outra forma, obtemos um anel de polinômios
em uma variável diferente.

::

    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<y> = PolynomialRing(QQ)
    sage: x == y
    False
    sage: R == S
    False
    sage: R(y)
    x
    sage: R(y^2 - 17)
    x^2 - 17

O anel é determinado pela variável. Note que criar um outro anel com
variável indeterminada ``x`` não retorna um anel diferente.

::

    sage: R = PolynomialRing(QQ, "x")
    sage: T = PolynomialRing(QQ, "x")
    sage: R == T
    True      
    sage: R is T
    True
    sage: R.0 == T.0
    True

O Sage também possui suporte para séries de potências e séries de
Laurent sobre um anel arbitrário. No seguinte exemplo, nós criamos um
elemento de :math:`\GF{7}[[T]]` e dividimos para criar um elemento de
:math:`\GF{7}((T))`.

::

    sage: R.<T> = PowerSeriesRing(GF(7)); R
    Power Series Ring in T over Finite Field of size 7
    sage: f = T  + 3*T^2 + T^3 + O(T^4)
    sage: f^3
    T^3 + 2*T^4 + 2*T^5 + O(T^6)
    sage: 1/f
    T^-1 + 4 + T + O(T^2)
    sage: parent(1/f)
    Laurent Series Ring in T over Finite Field of size 7

Você também pode criar anéis de polinômios usando a notação de
colchetes duplos:

::

    sage: GF(7)[['T']]
    Power Series Ring in T over Finite Field of size 7

Polinômios em Mais De Uma Variável
----------------------------------

Para trabalhar com polinômios em várias variáveis, nós primeiro
declaramos o anel de polinômios e as variáveis.

::

    sage: R = PolynomialRing(GF(5),3,"z") # here, 3 = number of variables
    sage: R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Da mesma forma como ocorre com polinômios em uma variável, existem
três maneiras de fazer isso:

::

    sage: GF(5)['z0, z1, z2']
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5
    sage: R.<z0,z1,z2> = GF(5)[]; R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Se você quiser usar os nomes das variáveis com apenas uma letra, então
você pode usar os seguinte comando:

::

    sage: PolynomialRing(GF(5), 3, 'xyz')
    Multivariate Polynomial Ring in x, y, z over Finite Field of size 5

A seguir fazemos um pouco de aritmética.

::

    sage: z = GF(5)['z0, z1, z2'].gens()
    sage: z
    (z0, z1, z2)
    sage: (z[0]+z[1]+z[2])^2
    z0^2 + 2*z0*z1 + z1^2 + 2*z0*z2 + 2*z1*z2 + z2^2

Você também pode usar uma notação mais matemática para criar um anel
de polinômios.

::

    sage: R = GF(5)['x,y,z']
    sage: x,y,z = R.gens()
    sage: QQ['x']
    Univariate Polynomial Ring in x over Rational Field
    sage: QQ['x,y'].gens()
    (x, y)
    sage: QQ['x'].objgens()
    (Univariate Polynomial Ring in x over Rational Field, (x,))

Polinômios em mais de uma variável são implementados no Sage usando
dicionários em Python e a "representação distribuída" de um polinômio.
O Sage usa o Singular [Si]_, por exemplo, para o cálculo do maior
divisor comum e bases de Gröbner para ideais algébricos.

::

    sage: R, (x, y) = PolynomialRing(RationalField(), 2, 'xy').objgens()
    sage: f = (x^3 + 2*y^2*x)^2
    sage: g = x^2*y^2
    sage: f.gcd(g)
    x^2

A seguir criamos o ideal :math:`(f,g)` gerado por :math:`f` e
:math:`g`, simplesmente multiplicando ``(f,g)`` por ``R`` (nós
poderíamos também escrever ``ideal([f,g])`` ou ``ideal(f,g)``).

.. link

::

    sage: I = (f, g)*R; I
    Ideal (x^6 + 4*x^4*y^2 + 4*x^2*y^4, x^2*y^2) of Multivariate Polynomial 
    Ring in x, y over Rational Field
    sage: B = I.groebner_basis(); B
    [x^6, x^2*y^2]
    sage: x^2 in I
    False

A base de Gröbner acima não é uma lista mas sim uma sequência
imutável. Isso implica que ela possui universo (universe) e parente
(parent), e não pode ser modificada (o que é bom pois ocasionaria
erros em outras rotinas que usam bases de Gröbner).

.. link

::

    sage: B.universe()
    Multivariate Polynomial Ring in x, y over Rational Field
    sage: B[1] = x
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Um pouco (não tanto quanto gostaríamos) de álgebra comutativa está
disponível no Sage, implementado via Singular. Por exemplo, podemos
calcular a decomposição primaria e primos associados de :math:`I`:

.. link

::

    sage: I.primary_decomposition()
    [Ideal (x^2) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y^2, x^6) of Multivariate Polynomial Ring in x, y over Rational Field]
    sage: I.associated_primes()
    [Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field]
