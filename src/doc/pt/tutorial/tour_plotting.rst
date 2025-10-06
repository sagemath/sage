.. sage-doctest: needs sage.plot sage.symbolic

.. _section-plot:

Gráficos
========

O Sage pode produzir gráficos bidimensionais e tridimensionais.

Gráficos Bidimensionais
-----------------------

Em duas dimensões, o Sage pode desenhar círculos, linhas, e polígonos;
gráficos de funções em coordenadas retangulares, e também coordenadas
polares; gráficos de contorno e gráficos de campos vetoriais.
Apresentamos alguns exemplos desses gráficos aqui. Para mais exemplos
de gráficos com o Sage, veja :ref:`section-systems` e
:ref:`section-maxima`, e também a documentação `Sage Constructions
<http://doc.sagemath.org/html/en/constructions/>`_.

Este comando produz um círculo amarelo de raio 1, centrado na origem.

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0))
    Graphics object consisting of 1 graphics primitive

Você pode também produzir um círculo preenchido:

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0), fill=True)
    Graphics object consisting of 1 graphics primitive

Outra possibilidade é criar um círculo atribuindo-o a uma variável;
isso não cria um gráfico:

::

    sage: c = circle((0,0), 1, rgbcolor=(1,1,0))

Para criar o gráfico, use ``c.show()`` ou ``show(c)``, da seguinte
forma:

.. link

::

    sage: c.show()

Alternativamente, o comando ``c.save('filename.png')`` salva o gráfico
no arquivo citado.

Agora, esses 'círculos' parecem mais elipses porque os eixos estão em
escalas diferentes. Você pode alterar isso:

.. link

::

    sage: c.show(aspect_ratio=1)

O comando ``show(c, aspect_ratio=1)`` produz o mesmo resultado, ou
você pode salvar a figura usando ``c.save('filename.png',
aspect_ratio=1)``.

É fácil criar o gráfico de funções simples:

::

    sage: plot(cos, (-5,5))
    Graphics object consisting of 1 graphics primitive

Após especificar uma variável, você também pode criar gráficos
paramétricos:

::

    sage: x = var('x')
    sage: parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    Graphics object consisting of 1 graphics primitive

É importante notar que os eixos dos gráficos vão se intersectar apenas
se a origem estiver no escopo do gráfico, e que valores grandes podem
ser representados usando notação científica.

::

    sage: plot(x^2,(x,300,500))
    Graphics object consisting of 1 graphics primitive

Você pode combinar vários gráficos somando-os:

::

    sage: x = var('x')
    sage: p1 = parametric_plot((cos(x),sin(x)),(x,0,2*pi),rgbcolor=hue(0.2))
    sage: p2 = parametric_plot((cos(x),sin(x)^2),(x,0,2*pi),rgbcolor=hue(0.4))
    sage: p3 = parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    sage: show(p1+p2+p3, axes=false)

Uma boa forma de produzir figuras preenchidas é criar uma lista de
pontos (``L`` no exemplo abaixo) e então usar o comando ``polygon``
para fazer o gráfico do polígono formado por esses pontos. Por
exemplo, aqui está um "deltoid" verde:

::

    sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),
    ....: 2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,3/4,1/2))
    sage: p
    Graphics object consisting of 1 graphics primitive

Digite ``show(p, axes=false)`` para visualizar isso sem os eixos.

Você pode adicionar texto ao gráfico:

::

    sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),
    ....: 6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))
    sage: t = text("hypotrochoid", (5,4), rgbcolor=(1,0,0))
    sage: show(p+t)

Professores de cálculo frequentemente desenham o seguinte gráfico na
lousa: não apenas um ramo do arco-seno, mas vários deles: isto é, o
gráfico de :math:`y=\sin(x)` para :math:`x` entre :math:`-2\pi` e
:math:`2\pi`, refletido com respeito a reta :math:`x=y`. Os seguintes
comandos fazem isso:

::

    sage: v = [(sin(x),x) for x in srange(-2*float(pi),2*float(pi),0.1)]
    sage: line(v)
    Graphics object consisting of 1 graphics primitive

Como a função tangente possui imagem maior do que o seno, se você usar
o mesmo método para fazer o gráfico da função inversa da função
tangente, você deve alterar as coordenadas mínima e máxima para o eixo
*x*:

::

    sage: v = [(tan(x),x) for x in srange(-2*float(pi),2*float(pi),0.01)]
    sage: show(line(v), xmin=-20, xmax=20)

O Sage também cria gráficos usando coordenadas polares, gráficos de
contorno e gráficos de campos vetoriais (para tipos especiais de
funções). Aqui está um exemplo de gráfico de contorno:

::

    sage: f = lambda x,y: cos(x*y)
    sage: contour_plot(f, (-4, 4), (-4, 4))
    Graphics object consisting of 1 graphics primitive

Gráficos Tridimensionais
------------------------

O Sage pode ser usado para criar gráficos tridimensionais. Tanto no
Sage Notebook, como no console (linha de comando), esses gráficos serão
exibidos usando o software de código aberto [ThreeJS]_, que permite girar
e ampliar a figura usando o mouse.

Use ``plot3d`` para criar o gráfico de uma função da forma `f(x, y) =
z`:

::

    sage: x, y = var('x,y')
    sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
    Graphics3d Object

Alternativamente, você pode usar ``parametric_plot3d`` para criar o
gráfico de uma superfície onde cada coordenada `x, y, z` é determinada
por uma função de uma ou duas variáveis (os parâmetros, tipicamente
`u` e `v`). O gráfico anterior pode ser representado parametricamente
na forma:

::

    sage: u, v = var('u, v')
    sage: f_x(u, v) = u
    sage: f_y(u, v) = v
    sage: f_z(u, v) = u^2 + v^2
    sage: parametric_plot3d([f_x, f_y, f_z], (u, -2, 2), (v, -2, 2))
    Graphics3d Object

A terceira forma de fazer um gráfico de uma superfície no Sage é
usando o comando ``implicit_plot3d``, que cria um gráfico de uma
superfície definida por uma equação `f(x, y, z) = 0` (isso define um
conjunto de pontos). Vamos fazer o gráfico de uma esfera usando a
expressão usual:

::

    sage: x, y, z = var('x, y, z')
    sage: implicit_plot3d(x^2 + y^2 + z^2 - 4, (x,-2, 2), (y,-2, 2), (z,-2, 2))
    Graphics3d Object

Aqui estão mais alguns exemplos:

`Yellow Whitney's umbrella <http://en.wikipedia.org/wiki/Whitney_umbrella>`__:

::

    sage: u, v = var('u,v')
    sage: fx = u*v
    sage: fy = u
    sage: fz = v^2
    sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -1, 1),
    ....: frame=False, color="yellow")
    Graphics3d Object

`Cross cap <http://en.wikipedia.org/wiki/Cross-cap>`__:

::

    sage: u, v = var('u,v')
    sage: fx = (1+cos(v))*cos(u)
    sage: fy = (1+cos(v))*sin(u)
    sage: fz = -tanh((2/3)*(u-pi))*sin(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....: frame=False, color="red")
    Graphics3d Object

Toro retorcido:

::

    sage: u, v = var('u,v')
    sage: fx = (3+sin(v)+cos(u))*cos(2*v)
    sage: fy = (3+sin(v)+cos(u))*sin(2*v)
    sage: fz = sin(u)+2*cos(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....: frame=False, color="red")
    Graphics3d Object

Lemniscata:

::

    sage: x, y, z = var('x,y,z')
    sage: f(x, y, z) = 4*x^2 * (x^2 + y^2 + z^2 + z) + y^2 * (y^2 + z^2 - 1)
    sage: implicit_plot3d(f, (x, -0.5, 0.5), (y, -1, 1), (z, -1, 1))
    Graphics3d Object
