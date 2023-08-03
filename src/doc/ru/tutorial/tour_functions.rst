.. _section-functions-issues:

Распространённые проблемы с функциями
=====================================

Некоторые аспекты определения функций (например, для дифференцирования
или построения графика) могут быть не ясны. В этом разделе мы обращаем
внимание на некоторые наиболее распространенные проблемы.

Далее показаны несколько способов определения того, что можно назвать
"функцией":

1. Определите функцию Python, как описано в разделе :ref:`section-functions`.
Для таких функций можно построить графики, но продифференцировать или
проинтегрировать их нельзя.

::

       sage: def f(z): return z^2
       sage: type(f)
       <... 'function'>
       sage: f(3)
       9
       sage: plot(f, 0, 2)
       Graphics object consisting of 1 graphics primitive

Обратите внимание на синтаксис в последней строчке. ``plot(f(z), 0, 2)``
выдаст ошибку, так как ``z`` - это переменная-болванка в определении ``f``,
которая не определена внутри данной конструкции. Просто ``f(z)`` возвратит
ошибку. Следующее будет работать в данном контексте, однако, в общем,
возникнут некоторые затруднения, но они могут быть проигнорированы (см. пункт 4).

.. link

::

       sage: var('z')   # определение переменной z для символьных вычислений
       z
       sage: f(z)
       z^2
       sage: plot(f(z), 0, 2)
       Graphics object consisting of 1 graphics primitive

В этом случае ``f(z)`` - это символьное выражение.

2. Определим "вызываемое символьное выражение". Оно может быть
продифференцировано, проинтегрировано, а также можно построить его график.

::

       sage: g(x) = x^2
       sage: g        # g отображает x в x^2
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

Если ``g`` — это вызываемое символьное выражение, ``g(x)`` — это
связянный с ним объект, но другого вида, для которого можно построить
график и который можно дифференциировать и т.д.

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

3. Можно использовать уже определенную функцию Sage — 'функцию исчисления'.
Для нее может быть построен график, она может быть продифференцирована
и проинтегрирована.

::

       sage: type(sin)
       <class 'sage.functions.trig.Function_sin'>
       sage: plot(sin, 0, 2)
       Graphics object consisting of 1 graphics primitive
       sage: type(sin(x))
       <class 'sage.symbolic.expression.Expression'>
       sage: plot(sin(x), 0, 2)
       Graphics object consisting of 1 graphics primitive

Сама по себе функция ``sin`` не может быть продифференцирована, по крайней
мере, не может произвести ``cos``.

::

       sage: f = sin
       sage: f.derivative()
       Traceback (most recent call last):
       ...
       AttributeError: ...

Использование ``f = sin(x)`` вместо ``sin`` работает, но лучше использовать
``f(x) = sin(x)`` для того, чтобы определить вызываемое символьное выражение.

::

       sage: S(x) = sin(x)
       sage: S.derivative()
       x |--> cos(x)

Далее следуют некоторые общие проблемы с объяснением:

\4. Случайная оценка.

::

       sage: def h(x):
       ....:     if x<2:
       ....:         return 0
       ....:     else:
       ....:         return x-2

Проблема: ``plot(h(x), 0, 4)`` построит кривую `y=x-2`.
Причина: В команде ``plot(h(x), 0, 4)`` сначала оценивается ``h(x)``,
что означает подставку ``x`` в функцию ``h`` и оценку ``x<2``.

.. link

::

       sage: type(x<2)
       <class 'sage.symbolic.expression.Expression'>

Решение: Не используйте ``plot(h(x), 0, 4)``; используйте:

.. link

::

       sage: plot(h, 0, 4)
       Graphics object consisting of 1 graphics primitive

\5. Ошибочное создание константы вместо функции.

::

       sage: f = x
       sage: g = f.derivative()
       sage: g
       1

Проблема: ``g(3)``, например, возвратит ошибку с сообщением
"ValueError: the number of arguments must be less than or equal to 0."

.. link

::

       sage: type(f)
       <class 'sage.symbolic.expression.Expression'>
       sage: type(g)
       <class 'sage.symbolic.expression.Expression'>

``g`` не является функцией, это константа, поэтому она не имеет
переменных, и вы можете вставлять что угодно в нее.

Решение: есть несколько возможных путей.

- Определить ``f`` изначально как символьное выражение.

::

         sage: f(x) = x        # вместо 'f = x'
         sage: g = f.derivative()
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>

- Либо вместе с ``f``, определенной выше, определить ``g`` как символьное выражение.

::

         sage: f = x
         sage: g(x) = f.derivative()  # вместо 'g = f.derivative()'
         sage: g
         x |--> 1
         sage: g(3)
         1
         sage: type(g)
         <class 'sage.symbolic.expression.Expression'>

- Либо с ``f`` и ``g``, заданными, как показано выше, создать переменную,
  под которую подставляются значения.

::

         sage: f = x
         sage: g = f.derivative()
         sage: g
         1
         sage: g(x=3)    # вместо 'g(3)'
         1

Есть еще один способ, как определить различие между производными
``f = x`` и ``f(x) = x``

::

       sage: f(x) = x
       sage: g = f.derivative()
       sage: g.variables()  # переменные, которые присутствуют в g
       ()
       sage: g.arguments()  # аргументы, которые могут быть подставлены в g
       (x,)
       sage: f = x
       sage: h = f.derivative()
       sage: h.variables()
       ()
       sage: h.arguments()
       ()

Как показывает данный пример, ``h`` не принимает аргументов,
поэтому ``h(3)`` вернет ошибку.
