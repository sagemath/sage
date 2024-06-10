.. _section-functions-issues:

Algumas Questões Frequentes sobre Funções
=========================================

Alguns aspectos sobre definição de funções (por exemplo, para
diferenciação, ou para criar gráficos) podem se tornar confusos. Nesta
seção, procuramos tratar algumas questões relevantes.

Aqui estão várias formas de definir objetos que merecem ser chamamos
de "funções":

1. Defina uma função em Python, como descrito em
:ref:`section-functions`. Essas funções podem ser usadas para criar
gráficos, mas não podem ser diferenciadas ou integradas.

::

       sage: def f(z): return z^2
       sage: type(f)
       <... 'function'>
       sage: f(3)
       9
       sage: plot(f, 0, 2)
       Graphics object consisting of 1 graphics primitive

Na última linha, observe a sintaxe. Se fosse usado ``plot(f(z), 0,
2)`` ocorreria um erro, porque ``z`` é uma variável muda na definição
de ``f`` e não está definida fora do contexto da função. De fato,
somente ``f(z)`` já provoca um erro. Os seguintes comandos vão
funcionar neste caso, embora em geral eles devam ser evitados pois
podem ocasionar erros (veja o item 4 abaixo).

.. link

::

       sage: var('z')   # define z to be a variable
       z
       sage: f(z)
       z^2
       sage: plot(f(z), 0, 2)
       Graphics object consisting of 1 graphics primitive

Acima, ``f(z)`` é uma expressão simbólica, o próximo item na nossa
lista.

2. Defina um "expressão simbólica que pode ser evocada". Essas podem
ser usadas para criar gráficos, e podem ser diferenciadas ou
integradas.

::

       sage: g(x) = x^2
       sage: g        # g sends x to x^2
       x |--> x^2
       sage: g(3)
       9
       sage: Dg = g.derivative(); Dg
       x |--> 2*x
       sage: Dg(3)
       6
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(g, 0, 2)
       Graphics object consisting of 1 graphics primitive

Note que enquanto ``g`` é uma expressão simbólica que pode ser
evocada, ``g(x)`` é um objeto diferente, embora relacionado, que pode
ser usado para criar gráficos, ou ser diferenciado, integrado, etc.,
embora com algumas ressalvas: veja o item 5 abaixo.

.. link

::

       sage: g(x)
       x^2
       sage: type(g(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: g(x).derivative()
       2*x
       sage: plot(g(x), 0, 2)
       Graphics object consisting of 1 graphics primitive

3. Use uma função pré-definida. Essas podem ser representadas em
gráfico, e com uma pequena ajuda, diferenciadas e integradas.

::

       sage: type(sin)
       <class 'sage.functions.trig.Function_sin'>
       sage: plot(sin, 0, 2)
       Graphics object consisting of 1 graphics primitive
       sage: type(sin(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(sin(x), 0, 2)
       Graphics object consisting of 1 graphics primitive
       
Por si só, ``sin`` não pode ser diferenciado, pelo menos não para
produzir ``cos``.

::

       sage: f = sin
       sage: f.derivative()
       Traceback (most recent call last):
       ...
       AttributeError: ...

Usando ``f = sin(x)`` no lugar de ``sin`` funciona, mas é ainda melhor
usar ``f(x) = sin(x)`` para definir uma expressão simbólica que pode
ser evocada.

::
   
       sage: S(x) = sin(x)
       sage: S.derivative()
       x |--> cos(x)
       
Aqui estão alguns problemas comuns, com explicações:

\4. Cálculo acidental.

::

       sage: def h(x):
       ....:     if x<2:
       ....:         return 0
       ....:     else:
       ....:         return x-2

O problema: ``plot(h(x), 0, 4)`` cria o gráfico da reta `y=x-2`, não
da função definida por ``h``. O motivo? No comando ``plot(h(x), 0,
4)``, primeiro ``h(x)`` é calculada: isso significa substituir ``x``
na função ``h``, o que significa que ``x<2`` é calculado.

.. link

::

       sage: type(x<2)
       <class 'sage.symbolic.expression.Expression'>

Quando uma equação simbólica é calculada, como na definição de ``h``,
se ela não é obviamente verdadeira, então ela retorna False. Logo
``h(x)`` é calculada como ``x-2``, e essa é a função que será
representada no gráfico.

A solução: não use ``plot(h(x), 0, 4)``; em vez disso, use

.. link

::

       sage: plot(h, 0, 4)
       Graphics object consisting of 1 graphics primitive

\5. Acidentalmente produzindo uma constante em vez de uma função.

::

       sage: f = x
       sage: g = f.derivative() 
       sage: g
       1

O problema: ``g(3)``, por exemplo, retorna o erro "ValueError: the
number of arguments must be less than or equal to 0."

.. link

::

       sage: type(f)
       <class 'sage.symbolic.expression.Expression'>
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>
       
``g`` não é uma função, é uma constante, logo não possui variáveis
associadas, e você não pode substituir nenhum valor em ``g``.

Solução: existem vária opções.

- Defina ``f`` inicialmente como uma expressão simbólica.

::

         sage: f(x) = x        # instead of 'f = x'
         sage: g = f.derivative()
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>

- Ou com ``f`` como definida originalmente, defina ``g`` como uma
  expressão simbólica.

::

         sage: f = x
         sage: g(x) = f.derivative()  # instead of 'g = f.derivative()'
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>

- Ou com ``f`` e ``g`` como definidas originalmente, especifique a
  variável para a qual você está substituindo.

::

         sage: f = x
         sage: g = f.derivative()
         sage: g
         1
         sage: g(x=3)    # instead of 'g(3)'
         1

Finalmente, aqui vai mais uma forma de saber a diferença entre as
derivadas de ``f = x`` e ``f(x) = x``.

::

       sage: f(x) = x 
       sage: g = f.derivative()
       sage: g.variables()  # the variables present in g
       ()
       sage: g.arguments()  # the arguments which can be plugged into g
       (x,)
       sage: f = x
       sage: h = f.derivative()
       sage: h.variables()
       ()
       sage: h.arguments()
       ()
       
Como esse exemplo procura ilustrar, ``h`` não aceita argumentos, e é
por isso que ``h(3)`` retorna um erro.
