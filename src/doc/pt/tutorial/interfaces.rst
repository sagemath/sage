.. linkall

**********
Interfaces
**********

Uma característica central do Sage é que ele permite fazer cálculos
com objetos em vários sistemas de álgebra computacional usando uma
interface comum e uma linguagem de programação clara.

Os métodos console e interact de uma interface executam tarefas bem
diferentes. Por exemplo, usando GAP:

#. ``gap.console()``: Isso abre um console do GAP - o controle é
   transferido para o GAP. Aqui o Sage não é nada mais do que uma
   forma conveniente de executar um programa, similar à shell Bash do
   Linux.

#. ``gap.interact()``: Essa é uma forma de interagir com uma instância
   do GAP que pode estar cheia de objetos do Sage. Você pode
   importar objetos nessa seção do GAP (até mesmo a partir da
   interface interativa), etc.


.. index: PARI; GP

GP/PARI
=======

O PARI é um programa em C muito compacto, maduro, e extremamente
otimizado cujo foco primário é teoria de números. Existem duas
interfaces distintas que podem ser usadas no Sage:


-  ``gp`` -- o "**G** do **P** ARI" interpretador, e

-  ``pari`` -- a biblioteca C do PARI.

Por exemplo, os seguintes comandos são duas formas de realizar a mesma
coisa. Eles parecem idênticos, mas o resultado é na verdade
diferente, e o que acontece por trás da cena é bastante diferente.

::

    sage: gp('znprimroot(10007)')
    Mod(5, 10007)
    sage: pari('znprimroot(10007)')
    Mod(5, 10007)

No primeiro caso, uma cópia separada do interpretador GP é iniciada
como um servidor, e a string ``znprimroot(10007)`` é enviada,
calculada pelo GP, e o resultado é armazenado em uma variável no GP
(que ocupa espaço na memória dos processos do GP que não serão
liberados). Então o valor dessa variável é exibido. No segundo caso,
nenhum programa separado é iniciado, e a string ``znprimroot(10007)``
é calculada por uma certa função da biblioteca C do PARI. O resultado
é armazenado na memória em uso pelo Python, que é liberada quando a
variável não for mais referenciada. Os objetos possuem tipos
diferentes:

::

    sage: type(gp('znprimroot(10007)'))
    <class 'sage.interfaces.gp.GpElement'>
    sage: type(pari('znprimroot(10007)'))
    <class 'cypari2.gen.Gen'>

Então qual eu devo usar? Depende do que você está fazendo. A interface
GP pode fazer absolutamente tudo o que você poderia fazer na linha de
comando do GP/PARI, pois está simplesmente executando esse programa.
Em particular, você pode carregar programas complicados em PARI e
executá-los. Por outro lado, a interface do PARI (via a biblioteca C)
é muito mais restritiva. Primeiro, nem todas as funções foram
implementadas. Segundo, bastante código, por exemplo, envolvendo
integração numérica, não irá funcionar através da interface PARI.
Todavia, a interface PARI pode ser significamente mais rápida e mais
robusta do que a interface GP.

(Se a interface GP ficar sem memória para calcular algum comando, ela
irá silenciosamente e automaticamente duplicar a memória alocada e
repetir o comando solicitado. Então os seus cálculos não irão ser
interrompidos se você não antecipou corretamente a quantidade de
memória que seria necessária. Esse é um truque útil que a interface
usual do GP não parece fornecer. Com respeito à interface da
biblioteca C do PARI, ela imediatamente copia cada objeto criado para
fora da pilha de memória, logo essa pilha nunca irá crescer. Contudo,
cada objeto não deve exceder 100MB de tamanho, ou a pilha irá estourar
quando o objeto for criado. Essa procedimento de cópia impõe uma leve
pena sobre a performace.)

Em resumo, o Sage usa a biblioteca C do pari para fornecer
funcionalidade similar à fornecida pelo interpretador GP/PARI, exceto
que com um gerenciamento sofisticado de memória e a linguagem de
programação Python.

Primeiro criamos uma lista do PARI a partir de uma lista do Python.

::

    sage: v = pari([1,2,3,4,5])
    sage: v
    [1, 2, 3, 4, 5]
    sage: type(v)
    <class 'cypari2.gen.Gen'>

Cada objeto do PARI é do tipo ``Gen``. O tipo PARI do objeto
subjacente pode ser obtido usando a função ``type``.

::

    sage: v.type()
    't_VEC'

Em PARI, para criar uma curva elíptica digitamos
``ellinit([1,2,3,4,5])``. Em Sage é similar, exceto que ``ellinit`` é
um método que pode ser chamado em qualquer objeto do PARI, por
exemplo, ``t_VEC`` `v`.

::

    sage: e = v.ellinit()
    sage: e.type()
    't_VEC'
    sage: pari(e)[:13]
    [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351]

Agora que temos um objeto de curva elíptica, podemos calcular algumas
coisas a respeito dele.

::

    sage: e.elltors()
    [1, [], []]
    sage: e.ellglobalred()
    [10351, [1, -1, 0, -1], 1, [11, 1; 941, 1], [[1, 5, 0, 1], [1, 5, 0, 1]]]
    sage: f = e.ellchangecurve([1,-1,0,-1])
    sage: f[:5]
    [1, -1, 0, 4, 3]

.. index: GAP

.. _section-gap:

GAP
===

O Sage vem com o GAP para matemática discreta computacional,
especialmente teoria de grupos.

Aqui está um exemplo com a função ``IdGroup`` do GAP, a qual usa a
base de dados opcional sobre grupos que precisa ser instalada
separadamente, como explicado abaixo.

::

    sage: G = gap('Group((1,2,3)(4,5), (3,4))')
    sage: G
    Group( [ (1,2,3)(4,5), (3,4) ] )
    sage: G.Center()
    Group( () )
    sage: G.IdGroup()
    [ 120, 34 ]
    sage: G.Order()
    120

Podemos realizar os mesmos cálculos no Sage sem explicitamente evocar
a interface do GAP da seguinte forma:

::

    sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
    sage: G.center()
    Subgroup generated by [()] of (Permutation Group with generators [(3,4), (1,2,3)(4,5)])
    sage: G.group_id()
    [120, 34]
    sage: n = G.order(); n
    120

Para algumas funcionalidades do GAP, deve-se instalar um pacote
Sage opcional. Isso pode ser feito com o comando::

    sage -i gap_packages


Singular
========

O Singular fornece uma biblioteca massiva e madura para bases de
Gröbner, máximo divisor comum para polinômios em várias variáveis,
bases de espaços de Riemann-Roch de uma curva plana, e fatorização,
entre outras coisas. Vamos ilustrar a fatorização de polinômios em
várias variáveis usando a interface do Sage para o Singular (não
digite ``...``):

::

    sage: R1 = singular.ring(0, '(x,y)', 'dp')
    sage: R1
    polynomial ring, over a field, global ordering
    //   coefficients: QQ
    //   number of vars : 2
    //        block   1 : ordering dp
    //                  : names    x y
    //        block   2 : ordering C
    sage: f = singular('9*y^8 - 9*x^2*y^7 - 18*x^3*y^6 - 18*x^5*y^6 +'
    ....:    '9*x^6*y^4 + 18*x^7*y^5 + 36*x^8*y^4 + 9*x^10*y^4 - 18*x^11*y^2 -'
    ....:    '9*x^12*y^3 - 18*x^13*y^2 + 9*x^16')

Agora que definimos :math:`f`, vamos imprimi-lo e fatorá-lo.

::

    sage: f
    9*x^16-18*x^13*y^2-9*x^12*y^3+9*x^10*y^4-18*x^11*y^2+36*x^8*y^4+18*x^7*y^5-18*x^5*y^6+9*x^6*y^4-18*x^3*y^6-9*x^2*y^7+9*y^8
    sage: f.parent()
    Singular
    sage: F = f.factorize(); F
    [1]:
       _[1]=9
       _[2]=x^6-2*x^3*y^2-x^2*y^3+y^4
       _[3]=-x^5+y^2
    [2]:
       1,1,2
    sage: F[1][2]
    x^6-2*x^3*y^2-x^2*y^3+y^4

Como com o exemplo para o GAP em :ref:`section-gap`, podemos calcular
a fatorização acima sem explicitamente usar a inteface do Singular
(todavia, implicitamente o Sage usa a interface do Singular para os
cálculos). Não digite ``...``:

::

    sage: x, y = QQ['x, y'].gens()
    sage: f = 9*y^8 - 9*x^2*y^7 - 18*x^3*y^6 - 18*x^5*y^6 + 9*x^6*y^4 \
    ....:     + 18*x^7*y^5 + 36*x^8*y^4 + 9*x^10*y^4 - 18*x^11*y^2 - 9*x^12*y^3 \
    ....:     - 18*x^13*y^2 + 9*x^16
    sage: factor(f)
    (9) * (-x^5 + y^2)^2 * (x^6 - 2*x^3*y^2 - x^2*y^3 + y^4)

.. _section-maxima:

Maxima
======

O Maxima está incluido no Sage, assim como uma implementação do Lisp.
O pacote gnuplot (que o Maxima usa para criar gráficos) é distribuído
como um pacote adicional do Sage. Entre outras coisas, o Maxima
executa manipulações simbólicas. Ele pode integrar e diferenciar
funções simbolicamente, resolver EDOs de primeira ordem, grande parte
das EDOs lineares de segunda ordem, e tem implementado o método da
transformada de Laplace para EDOs lineares de qualquer ordem. O Maxima
também suporta uma série de funções especiais, é capaz de criar
gráficos via gnuplot, e possui métodos para resolver equações
polinômiais e manipular matrizes (por exemplo, escalonar e calcular
autovalores e autovetores).

Nós ilustramos a interface Sage/Maxima construíndo uma matriz cuja
entrada :math:`i,j` é :math:`i/j`, para :math:`i,j=1,\ldots,4`.

::

    sage: f = maxima.eval('ij_entry[i,j] := i/j')
    sage: A = maxima('genmatrix(ij_entry,4,4)'); A
    matrix([1,1/2,1/3,1/4],[2,1,2/3,1/2],[3,3/2,1,3/4],[4,2,4/3,1])
    sage: A.determinant()
    0
    sage: A.echelon()
    matrix([1,1/2,1/3,1/4],[0,0,0,0],[0,0,0,0],[0,0,0,0])
    sage: A.eigenvalues()
    [[0,4],[3,1]]
    sage: A.eigenvectors().sage()
    [[[0, 4], [3, 1]], [[[1, 0, 0, -4], [0, 1, 0, -2], [0, 0, 1, -4/3]], [[1, 2, 3, 4]]]]

Aqui vai outro exemplo:

::

    sage: A = maxima("matrix ([1, 0, 0], [1, -1, 0], [1, 3, -2])")
    sage: eigA = A.eigenvectors()
    sage: V = VectorSpace(QQ,3)
    sage: eigA
    [[[-2,-1,1],[1,1,1]],[[[0,0,1]],[[0,1,3]],[[1,1/2,5/6]]]]
    sage: v1 = V(sage_eval(repr(eigA[1][0][0]))); lambda1 = eigA[0][0][0]
    sage: v2 = V(sage_eval(repr(eigA[1][1][0]))); lambda2 = eigA[0][0][1]
    sage: v3 = V(sage_eval(repr(eigA[1][2][0]))); lambda3 = eigA[0][0][2]

    sage: M = MatrixSpace(QQ,3,3)
    sage: AA = M([[1,0,0],[1, - 1,0],[1,3, - 2]])
    sage: b1 = v1.base_ring()
    sage: AA*v1 == b1(lambda1)*v1
    True
    sage: b2 = v2.base_ring()
    sage: AA*v2 == b2(lambda2)*v2
    True
    sage: b3 = v3.base_ring()
    sage: AA*v3 == b3(lambda3)*v3
    True

Por fim, apresentamos um exemplo de como usar o Sage para criar
gráficos usando ``openmath``. Alguns desses exemplos são modificações
de exemplos do manual de referência do Maxima.

Um gráfico em duas dimensões de diversas funções (não digite ``...``):

::

    sage: maxima.plot2d('[cos(7*x),cos(23*x)^4,sin(13*x)^3]', '[x,0,1]',  # not tested
    ....:               '[plot_format,openmath]')

Um gráfico em 3D que você pode mover com o seu mouse:

::

    sage: maxima.plot3d("2^(-u^2 + v^2)", "[u, -3, 3]", "[v, -2, 2]",  # not tested
    ....:               '[plot_format, openmath]')

    sage: maxima.plot3d("atan(-x^2 + y^3/4)", "[x, -4, 4]", "[y, -4, 4]",  # not tested
    ....:               "[grid, 50, 50]",'[plot_format, openmath]')

O próximo gráfico é a famosa faixa de Möbious:

::

    sage: maxima.plot3d("[cos(x)*(3 + y*cos(x/2)), sin(x)*(3 + y*cos(x/2)),"  # not tested
    ....:               "y*sin(x/2)]", "[x, -4, 4]", "[y, -4, 4]",
    ....:               '[plot_format, openmath]')

E agora a famosa garrafa de Klein:

::

    sage: _ = maxima("expr_1: 5*cos(x)*(cos(x/2)*cos(y) + sin(x/2)*sin(2*y)+ 3.0)"
    ....:        "- 10.0")
    sage: _ = maxima("expr_2: -5*sin(x)*(cos(x/2)*cos(y) + sin(x/2)*sin(2*y)+ 3.0)")
    sage: _ = maxima("expr_3: 5*(-sin(x/2)*cos(y) + cos(x/2)*sin(2*y))")
    sage: maxima.plot3d("[expr_1, expr_2, expr_3]", "[x, -%pi, %pi]", # not tested
    ....:               "[y, -%pi, %pi]", "['grid, 40, 40]",
    ....:               '[plot_format, openmath]')
