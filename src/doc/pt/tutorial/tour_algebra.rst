Álgebra Elementar e Cálculo
===========================

O Sage pode realizar diversos cálculos em álgebra elementar e cálculo
diferencial e integral: por exemplo, encontrar soluções de equações,
diferenciar, integrar, e calcular a transformada de Laplace. Veja a
documentação em `Sage Constructions
<http://doc.sagemath.org/html/en/constructions/>`_ para mais exemplos.

Resolvendo equações
-------------------

Resolvendo equações exatamente
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A função ``solve`` resolve equações. Para usá-la, primeiro especifique
algumas variáveis; então os argumentos de ``solve`` são uma equação
(ou um sistema de equações), juntamente com as variáveis para as
quais resolver:

::

    sage: x = var('x')
    sage: solve(x^2 + 3*x + 2, x)
    [x == -2, x == -1]

Você pode resolver equações para uma variável em termos das outras:

::

    sage: x, b, c = var('x b c')
    sage: solve([x^2 + b*x + c == 0],x)
    [x == -1/2*b - 1/2*sqrt(b^2 - 4*c), x == -1/2*b + 1/2*sqrt(b^2 - 4*c)]

Você pode resolver para diversas variáveis:

::

    sage: x, y = var('x, y')
    sage: solve([x+y==6, x-y==4], x, y)
    [[x == 5, y == 1]]

O seguinte exemplo, que mostra como usar o Sage para resolver um
sistema de equações não-lineares, foi sugerido por Jason Grout:
primeiro, resolvemos o sistemas simbolicamente:

::

    sage: var('x y p q')
    (x, y, p, q)
    sage: eq1 = p+q==9
    sage: eq2 = q*y+p*x==-6
    sage: eq3 = q*y^2+p*x^2==24
    sage: solve([eq1,eq2,eq3,p==1],p,q,x,y)
    [[p == 1, q == 8, x == -4/3*sqrt(10) - 2/3, y == 1/6*sqrt(10) - 2/3], [p == 1, q == 8, x == 4/3*sqrt(10) - 2/3, y == -1/6*sqrt(10) - 2/3]]

Para obter soluções numéricas aproximadas, podemos usar:

.. link

::

    sage: solns = solve([eq1,eq2,eq3,p==1],p,q,x,y, solution_dict=True)
    sage: [[s[p].n(30), s[q].n(30), s[x].n(30), s[y].n(30)] for s in solns]
    [[1.0000000, 8.0000000, -4.8830369, -0.13962039],
     [1.0000000, 8.0000000, 3.5497035, -1.1937129]]

(A função ``n`` imprime uma aproximação numérica, e o argumento é o
número de bits de precisão.)

Resolvendo Equações Numericamente
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frequentemente, ``solve`` não será capaz de encontrar uma solução
exata para uma equação ou sistema de equações. Nesse caso, você pode
usar ``find_root`` para encontrar uma solução numérica. Por exemplo,
``solve`` não encontra uma solução para a equação abaixo::

    sage: theta = var('theta')
    sage: solve(cos(theta)==sin(theta), theta)
    [sin(theta) == cos(theta)]

Por outro lado, podemos usar ``find_root`` para encontrar uma solução
para a equação acima no intervalo :math:`0 < \phi < \pi/2`::

    sage: phi = var('phi')
    sage: find_root(cos(phi)==sin(phi),0,pi/2)
    0.785398163397448...

Diferenciação, Integração, etc.
-------------------------------

O Sage é capaz de diferenciar e integrar diversas funções. Por
exemplo, para diferenciar :math:`\sin(u)` com respeito a :math:`u`,
faça o seguinte:

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

Para calcular a quarta derivada de :math:`\sin(x^2)`:

::

    sage: diff(sin(x^2), x, 4)
    16*x^4*sin(x^2) - 48*x^2*cos(x^2) - 12*sin(x^2)

Para calcular as derivadas parciais de :math:`x^2+17y^2` com respeito
a *x* e *y*, respectivamente:

::

    sage: x, y = var('x,y')
    sage: f = x^2 + 17*y^2
    sage: f.diff(x)
    2*x
    sage: f.diff(y)
    34*y

Passamos agora para integrais, tanto indefinidas como definidas. Para
calcular :math:`\int x\sin(x^2)\, dx` e :math:`\int_0^1
\frac{x}{x^2+1}\, dx`:

::

    sage: integral(x*sin(x^2), x)
    -1/2*cos(x^2)
    sage: integral(x/(x^2+1), x, 0, 1)
    1/2*log(2)

Para calcular a decomposição em frações parciais de
:math:`\frac{1}{x^2-1}`:

::

    sage: f = 1/((1+x)*(x-1))
    sage: f.partial_fraction(x)
    -1/2/(x + 1) + 1/2/(x - 1)

.. _section-systems:

Resolvendo Equações Diferenciais
--------------------------------

Você pode usar o Sage para investigar equações diferenciais
ordinárias. Para resolver a equação :math:`x'+x-1=0`:

::

    sage: t = var('t')    # define a variable t
    sage: x = function('x')(t)   # define x to be a function of that variable
    sage: DE = diff(x, t) + x - 1
    sage: desolve(DE, [x,t])
    (_C + e^t)*e^(-t)

Esse método usa a interface do Sage para o Maxima [Max]_. Logo, o
formato dos resultados é um pouco diferente de outros cálculos
realizados no Sage. Nesse caso, o resultado diz que a solução geral da
equação diferencial é :math:`x(t) = e^{-t}(e^{t}+c)`.

Você pode calcular a transformada de Laplace também; a transformada de
Laplace de :math:`t^2e^t -\sin(t)` é calculada da seguinte forma:

::

    sage: s = var("s")
    sage: t = var("t")
    sage: f = t^2*exp(t) - sin(t)
    sage: f.laplace(t,s)
    -1/(s^2 + 1) + 2/(s - 1)^3

A seguir, um exemplo mais complicado. O deslocamento, com respeito à
posição de equilíbrio, de duas massas presas a uma parede através de
molas, conforme a figura abaixo,

::

    |------\/\/\/\/\---|massa1|----\/\/\/\/\/----|massa2|
             mola1                    mola2

é modelado pelo sistema de equações diferenciais de segunda ordem

.. math::

    m_1 x_1'' + (k_1+k_2) x_1 - k_2 x_2 = 0

    m_2 x_2''+ k_2 (x_2-x_1) = 0,


onde, para :math:`i=1,2`, :math:`m_{i}` é a massa do objeto *i*,
:math:`x_{i}` é o deslocamento com respeito à posição de equilíbrio da
massa *i*, e :math:`k_{i}` é a constante de mola para a mola *i*.

**Exemplo:** Use o Sage para resolver o problema acima com
:math:`m_{1}=2`, :math:`m_{2}=1`, :math:`k_{1}=4`,
:math:`k_{2}=2`, :math:`x_{1}(0)=3`, :math:`x_{1}'(0)=0`,
:math:`x_{2}(0)=3`, :math:`x_{2}'(0)=0`.

Solução: Primeiramente, calcule a transformada de Laplace da primeira
equação (usando a notação :math:`x=x_{1}`, :math:`y=x_{2}`):

::

    sage: de1 = maxima("2*diff(x(t),t, 2) + 6*x(t) - 2*y(t)")
    sage: lde1 = de1.laplace("t","s"); lde1.sage()
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

O resultado é um pouco difícil de ler, mas diz que

.. math:: -2x'(0) + 2s^2*X(s) - 2sx(0) - 2Y(s) + 6X(s) = 0


(onde a transformada de Laplace de uma função em letra minúscula
:math:`x(t)` é a função em letra maiúscula :math:`X(s)`). Agora,
calcule a transformada de Laplace da segunda equação:

::

    sage: t,s = SR.var('t,s')
    sage: x = function('x')
    sage: y = function('y')
    sage: f = 2*x(t).diff(t,2) + 6*x(t) - 2*y(t)
    sage: f.laplace(t,s)
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

O resultado significa que

.. math:: -Y'(0) + s^2Y(s) + 2Y(s) - 2X(s) - sy(0) = 0.


Em seguida, substitua a condição inicial para :math:`x(0)`,
:math:`x'(0)`, :math:`y(0)`, e :math:`y'(0)`, e resolva as equações
resultantes:

::

    sage: var('s X Y')
    (s, X, Y)
    sage: eqns = [(2*s^2+6)*X-2*Y == 6*s, -2*X +(s^2+2)*Y == 3*s]
    sage: solve(eqns, X,Y)
    [[X == 3*(s^3 + 3*s)/(s^4 + 5*s^2 + 4),
      Y == 3*(s^3 + 5*s)/(s^4 + 5*s^2 + 4)]]

Agora calcule a transformada de Laplace inversa para obter a resposta:

::

    sage: var('s t')
    (s, t)
    sage: inverse_laplace((3*s^3 + 9*s)/(s^4 + 5*s^2 + 4),s,t)
    cos(2*t) + 2*cos(t)
    sage: inverse_laplace((3*s^3 + 15*s)/(s^4 + 5*s^2 + 4),s,t)
    -cos(2*t) + 4*cos(t)

Portanto, a solução é

.. math:: x_1(t) = \cos(2t) + 2\cos(t), \quad x_2(t) = 4\cos(t) - \cos(2t).


Ela pode ser representada em um gráfico parametricamente usando os
comandos

::

    sage: t = var('t')
    sage: P = parametric_plot((cos(2*t) + 2*cos(t), 4*cos(t) - cos(2*t) ),
    ....: (t, 0, 2*pi), rgbcolor=hue(0.9))
    sage: show(P)

As componentes individuais podem ser representadas em gráfico usando

::

    sage: t = var('t')
    sage: p1 = plot(cos(2*t) + 2*cos(t), (t,0, 2*pi), rgbcolor=hue(0.3))
    sage: p2 = plot(4*cos(t) - cos(2*t), (t,0, 2*pi), rgbcolor=hue(0.6))
    sage: show(p1 + p2)

Leia mais sobre gráficos em :ref:`section-plot`. Veja a seção 5.5 de
[NagleEtAl2004]_ (em inglês) para mais informações sobre equações
diferenciais.


Método de Euler para Sistemas de Equações Diferenciais
------------------------------------------------------

No próximo exemplo, vamos ilustrar o método de Euler para EDOs de
primeira e segunda ordem. Primeiro, relembramos a ideia básica para
equações de primeira ordem. Dado um problema de valor inicial da forma

.. math::

    y'=f(x,y), \quad y(a)=c,

queremos encontrar o valor aproximado da solução em :math:`x=b` com
:math:`b>a`.

Da definição de derivada segue que

.. math::  y'(x) \approx \frac{y(x+h)-y(x)}{h},


onde :math:`h>0` é um número pequeno. Isso, juntamente com a equação
diferencial, implica que :math:`f(x,y(x))\approx
\frac{y(x+h)-y(x)}{h}`. Agora resolvemos para :math:`y(x+h)`:

.. math::   y(x+h) \approx y(x) + h*f(x,y(x)).


Se chamarmos :math:`h f(x,y(x))` de "termo de correção", :math:`y(x)`
de "valor antigo de *y*", e :math:`y(x+h)` de "novo valor de *y*",
então essa aproximação pode ser reescrita como

.. math::   y_{novo} \approx y_{antigo} + h*f(x,y_{antigo}).


Se dividirmos o intervalo de *a* até *b* em *n* partes, de modo que
:math:`h=\frac{b-a}{n}`, então podemos construir a seguinte tabela.

============== ==================   ================
:math:`x`      :math:`y`            :math:`hf(x,y)`
============== ==================   ================
:math:`a`      :math:`c`            :math:`hf(a,c)`
:math:`a+h`    :math:`c+hf(a,c)`    ...
:math:`a+2h`   ...
...
:math:`b=a+nh` ???                  ...
============== ==================   ================


O objetivo é completar os espaços em branco na tabela, em uma linha
por vez, até atingirmos ???, que é a aproximação para :math:`y(b)`
usando o método de Euler.

A ideia para sistemas de EDOs é semelhante.

**Exemplo:** Aproxime numericamente :math:`z(t)` em :math:`t=1` usando
4 passos do método de Euler, onde :math:`z''+tz'+z=0`, :math:`z(0)=1`,
:math:`z'(0)=0`.

Devemos reduzir a EDO de segunda ordem a um sistema de duas EDOs de
primeira ordem (usando :math:`x=z`, :math:`y=z'`) e aplicar o método
de Euler:

::

    sage: t,x,y = PolynomialRing(RealField(10),3,"txy").gens()
    sage: f = y; g = -x - y * t
    sage: eulers_method_2x2(f,g, 0, 1, 0, 1/4, 1)
          t                x            h*f(t,x,y)                y       h*g(t,x,y)
          0                1                  0.00                0           -0.25
        1/4              1.0                -0.062            -0.25           -0.23
        1/2             0.94                 -0.12            -0.48           -0.17
        3/4             0.82                 -0.16            -0.66          -0.081
          1             0.65                 -0.18            -0.74           0.022

Portanto, :math:`z(1)\approx 0.65`.

Podemos também representar em um gráfico os pontos :math:`(x,y)` para
obter uma figura da solução aproximada. A função
``eulers_method_2x2_plot`` fará isso; para usá-la, precisamos definir
funções *f* e *g* que recebam um argumento com três coordenadas (*t*,
*x*, *y*).

::

    sage: f = lambda z: z[2]        # f(t,x,y) = y
    sage: g = lambda z: -sin(z[1])  # g(t,x,y) = -sin(x)
    sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)

A esta altura, ``P`` armazena dois gráficos: ``P[0]``, o gráfico de
*x* versus *t*, e ``P[1]``, o gráfico de *y* versus *t*. Podemos
visualizar os dois gráficos da seguinte forma:

.. link

::

    sage: show(P[0] + P[1])

(Para mais sobre gráficos, veja :ref:`section-plot`.)

Funções Especiais
-----------------

Diversos polinômios ortogonais e funções especiais estão
implementadas, usando tanto o PARI [GP]_ como o Maxima [Max]_. Isso
está documentado nas seções apropriadas ("Orthogonal polynomials" and
"Special functions", respectivamente) do manual de referência do Sage
(em inglês).

::

    sage: x = polygen(QQ, 'x')
    sage: chebyshev_U(2,x)
    4*x^2 - 1
    sage: bessel_I(1,1).n(250)
    0.56515910399248502720769602760986330732889962162109200948029448947925564096
    sage: bessel_I(1,1).n()
    0.56515910399248...
    sage: bessel_I(2,1.1).n()  # last few digits are random
    0.16708949925104...

No momento, essas funções estão disponíveis na interface do Sage
apenas para uso numérico. Para uso simbólico, use a interface do
Maxima diretamente, como no seguinte exemplo:

::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'
