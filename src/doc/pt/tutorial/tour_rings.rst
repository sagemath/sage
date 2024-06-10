.. _section-rings:

Anéis Básicos
=============

Quando se define matrizes, vetores, ou polinômios, é as vezes útil, e
as vezes necessário, especificar o "anel" sobre o qual o objeto será
definido. Um *anel* é uma estrutura matemática na qual se tem noções
de adição e multiplicação bem definidas; se você nunca ouviu falar
sobre anéis, você provavelmente só precisa saber a respeito dos
seguintes exemplos:

* os inteiros `\{..., -1, 0, 1, 2,... \}`, que são chamados ``ZZ`` no
  Sage.
* os números racionais -- i. e., frações, ou razões, de inteiros --
  que são chamados ``QQ`` no Sage.
* os números reais, chamados de ``RR`` no Sage.
* os números complexos, chamados de ``CC`` no Sage.

Você pode precisar saber sobre essas distinções porque o mesmo
polinômio, por exemplo, pode ser tratado diferentemente dependendo do
anel sobre o qual está definido. A propósito, o polinômio `x^2-2`
possui duas raízes, `\pm \sqrt{2}`. Essas raízes não são racionais,
logo, se você esta lidando com polinômios com coeficientes racionais,
os polinômios não serão fatorados. Com coeficientes reais, eles serão.
Portanto você pode querer especificar o anel para garantir que você
vai obter a informação que deseja. Os dois comandos a seguir definem
os conjuntos de polinômios com coeficientes racionais e coeficientes
reais, respectivamente. Os conjuntos são chamados "ratpoly" e
"realpoly", mas esses nomes não são importantes aqui; todavia, note
que as strings ".<t>" e ".<z>" especificam o nome das variáveis
usadas em cada caso.

::

    sage: ratpoly.<t> = PolynomialRing(QQ)
    sage: realpoly.<z> = PolynomialRing(RR)

Agora ilustramos a nossa discussão sobre fatorar `x^2-2`:

.. link

::

    sage: factor(t^2-2)
    t^2 - 2
    sage: factor(z^2-2)
    (z - 1.41421356237310) * (z + 1.41421356237310)

Comentários similares também se aplicam a matrizes: a forma reduzida
de uma matriz pode depender do anel sobre o qual ela esta definida,
como também pode os seus autovalores e autovetores. Para mais sobre
polinômios, veja :ref:`section-poly`, para mais sobre matrizes, veja
:ref:`section-linalg`.

O símbolo ``I`` representa a raiz quadrada de :math:`-1`; ``i`` é um
sinônimo de ``I``. Obviamente, isso não é um número racional::

    sage: i  # square root of -1
    I     
    sage: i in QQ
    False

Nota: O código acima pode não funcionar como esperado se a variável
``i`` estiver atribuída a um outro valor, por exemplo, se ela for
usada como a variável de um laço (loop). Nesse caso, digite::

    sage: reset('i')

para restabelecer o valor original de ``i``.

Há uma sutileza ao definir números complexos: como mencionado acima,
o símbolo ``i`` representa a raiz quadrada de `-1`, mas é uma raiz
quadrada de `-1` como número algébrico. Evocando ``CC(i)`` ou ``CC.0``
ou ``CC.gen(0)`` obtém-se a raiz de `-1` complexa. Aritmética
envolvendo tipos diferentes de números é possível graças ao que se
chama de coação, veja :ref:`section-coercion`.

::

    sage: i = CC(i)       # floating point complex number
    sage: i == CC.0
    True
    sage: a, b = 4/3, 2/3
    sage: z = a + b*i
    sage: z
    1.33333333333333 + 0.666666666666667*I
    sage: z.imag()        # imaginary part
    0.666666666666667
    sage: z.real() == a   # automatic coercion before comparison
    True
    sage: a + b
    2
    sage: 2*b == a
    True
    sage: parent(2/3)
    Rational Field
    sage: parent(4/2)
    Rational Field
    sage: 2/3 + 0.1       # automatic coercion before addition
    0.766666666666667
    sage: 0.1 + 2/3       # coercion rules are symmetric in Sage
    0.766666666666667

Aqui estão mais exemplos de anéis básicos em Sage. Como observado
acima, o anel dos números racionais pode ser referido usando ``QQ``,
ou também ``RationalField()`` (um *corpo*, ou *field* em inglês, é um
anel no qual a operação de multiplicação é comutativa, e todo elemento
não-nulo possui um elemento inverso com respeito à operação de
multiplicação. Logo, os racionais formam um corpo, mas os inteiros
não)::

    sage: RationalField()
    Rational Field
    sage: QQ
    Rational Field
    sage: 1/2 in QQ
    True

O número decimal ``1.2`` é considerado como um elemento de ``QQ``:
número decimais que são também racionais podem ser coagidos ao conjunto de
números racionais (veja :ref:`section-coercion`). Os números `\pi` e
`\sqrt{2}` não são racionais, todavia::

    sage: 1.2 in QQ
    True
    sage: pi in QQ
    False
    sage: pi in RR
    True
    sage: sqrt(2) in QQ
    False
    sage: sqrt(2) in CC
    True

Para uso em matemática mais avançada, o Sage também pode especificar
outros anéis, como corpos finitos, inteiros `p`-ádicos, o anel dos
números algébricos, anéis de polinômios, e anéis de matrizes. Aqui
está a construção de alguns deles::

    sage: GF(3)
    Finite Field of size 3
    sage: GF(27, 'a')  # need to name the generator if not a prime field
    Finite Field in a of size 3^3
    sage: Zp(5)
    5-adic Ring with capped relative precision 20
    sage: sqrt(3) in QQbar # algebraic closure of QQ
    True
