Álgebra Y Cálculo Básicos
==========================

Sage puede efectuar cómputos relacionados al algebra y cálculo básicos:
por ejemplo, encontrar soluciones de ecuaciones, diferenciación, integración y transformadas de Laplace.
Véa la documentación "Construcciones En Sage"  para más ejemplos.

Resolviendo Ecuaciones
----------------------

Resolviendo Ecuaciones De Manera Exacta
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

La función ``solve`` resuelve ecuaciones. Para usarla, primero no olvides especificar
algunas variables. Los argumentos de ``solve`` son una ecuación (o un
sistema de ecuaciones), junto con las variables a resolver:

::

    sage: x = var('x')
    sage: solve(x^2 + 3*x + 2, x)
    [x == -2, x == -1]

Puedes resolver ecuaciones en una variable respecto de las demás:

::

    sage: x, b, c = var('x b c')
    sage: solve([x^2 + b*x + c == 0],x)
    [x == -1/2*b - 1/2*sqrt(b^2 - 4*c), x == -1/2*b + 1/2*sqrt(b^2 - 4*c)]

Puedes también resolver ecuaciones en varias variables:

::

    sage: x, y = var('x, y')
    sage: solve([x+y==6, x-y==4], x, y)
    [[x == 5, y == 1]]

El siguiente ejemplo del uso de Sage para resolver un sistema de ecuaciones
no-lineales fue proporcionado por Jason Grout: primero, resolvemos el sistema
simbólicamente:

::

    sage: var('x y p q')
    (x, y, p, q)
    sage: eq1 = p+q==9
    sage: eq2 = q*y+p*x==-6
    sage: eq3 = q*y^2+p*x^2==24
    sage: solve([eq1,eq2,eq3,p==1],p,q,x,y)
    [[p == 1, q == 8, x == -4/3*sqrt(10) - 2/3, y == 1/6*sqrt(10) - 2/3], [p == 1, q == 8, x == 4/3*sqrt(10) - 2/3, y == -1/6*sqrt(10) - 2/3]]

Si queremos aproximaciones numéricas de las soluciones, podemos usar lo siguiente:

.. link

::

    sage: solns = solve([eq1,eq2,eq3,p==1],p,q,x,y, solution_dict=True)
    sage: [[s[p].n(30), s[q].n(30), s[x].n(30), s[y].n(30)] for s in solns]
    [[1.0000000, 8.0000000, -4.8830369, -0.13962039],
     [1.0000000, 8.0000000, 3.5497035, -1.1937129]]

(La función ``n`` imprime una aproximación numérica, y el
argumento es el número de bits de precisión.)

Resolviendo Ecuaciones Numéricamente
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A menudo, ``solve`` no podrá encontrar una solución exacta para
la ecuación o ecuaciones especificadas.  Cuando falla, puedes usar
``find_root`` para encontrar una solución numérica.  Por ejemplo, ``solve`` no
devuelve nada interesante para la siguiente ecuación::

    sage: theta = var('theta')
    sage: solve(cos(theta)==sin(theta), theta)
    [sin(theta) == cos(theta)]

Por otro lado, podemos usar ``find_root`` para encontrar una solución a la
ecuación de arriba en el rango :math:`0 < \theta < \pi/2`::

    sage: phi = var('phi')
    sage: find_root(cos(phi)==sin(phi),0,pi/2)
    0.785398163397448...

Diferenciación, Integración, etc.
----------------------------------

Sage sabe cómo diferenciar e integrar muchas funciones.
Por ejemplo, para diferenciar :math:`\sin(u)` con respecto a :math:`u`,
haz lo siguiente:

::

    sage: u = var('u')
    sage: diff(sin(u), u)
    cos(u)

Para calcular la cuarta derivada de :math:`\sin(x^2)`:

::

    sage: diff(sin(x^2), x, 4)
    16*x^4*sin(x^2) - 48*x^2*cos(x^2) - 12*sin(x^2)

Para calcular las derivadas parciales de :math:`x^2+17y^2` con
respecto a *x* e *y*, respectivamente:

::

    sage: x, y = var('x,y')
    sage: f = x^2 + 17*y^2
    sage: f.diff(x)
    2*x
    sage: f.diff(y)
    34*y

También podemos calcular integrales, tanto indefinidas como definidas.
Para calcular :math:`\int x\sin(x^2)\, dx` y :math:`\int_0^1 \frac{x}{x^2+1}\, dx`

::

    sage: integral(x*sin(x^2), x)
    -1/2*cos(x^2)
    sage: integral(x/(x^2+1), x, 0, 1)
    1/2*log(2)

Para calcular la descomposición en fracciones simples de
:math:`\frac{1}{x^2-1}`:

::

    sage: f = 1/((1+x)*(x-1))
    sage: f.partial_fraction(x)
    -1/2/(x + 1) + 1/2/(x - 1)

.. _section-systems:

Resolviendo Ecuaciones Diferenciales
------------------------------------

Puedes usar a Sage para investigar ecuaciones diferenciales ordinarias.
Para resolver la ecuación :math:`x'+x-1=0`:

::

    sage: t = var('t')    # defina una variable t
    sage: x = function('x')(t)   # defina x como una función de esa variable
    sage: DE = diff(x, t) + x - 1
    sage: desolve(DE, [x,t])
    (_C + e^t)*e^(-t)

Esto utiliza el interfaz a Maxima de Sage [Max]_, por lo que el resultado puede
diferir de otros resultados de Sage. En este caso, la salida nos dice que la
solución general a la ecuación diferencial es :math:`x(t) = e^{-t}(e^{t}+c)`.

También puedes calcular transformadas de Laplace; la transformada de Laplace
de :math:`t^2e^t -\sin(t)` se calcula como sigue:

::

    sage: s = var("s")
    sage: t = var("t")
    sage: f = t^2*exp(t) - sin(t)
    sage: f.laplace(t,s)
    -1/(s^2 + 1) + 2/(s - 1)^3

Veamos un ejemplo más complicado. El desplazamiento desde el punto de equilibrio
de dos resortes acoplados, sujetos a una pared a la izquierda

::

    |------\/\/\/\/\---|masa1|----\/\/\/\/\/----|masa2|
             resorte1               resorte2

está modelado por el sistema de ecuaciones diferenciales de segundo órden

.. math::

    m_1 x_1'' + (k_1+k_2) x_1 - k_2 x_2 = 0

    m_2 x_2''+ k_2 (x_2-x_1) = 0,

donde :math:`m_{i}` es la masa del objeto *i*, :math:`x_{i}` es
el desplazamiento desde el equilibrio de la masa *i*, y :math:`k_{i}`
es la constante de elasticidad del resorte *i*.

**Ejemplo:** Utiliza Sage para resolver el problema de arriba con
:math:`m_{1}=2`, :math:`m_{2}=1`, :math:`k_{1}=4`,
:math:`k_{2}=2`, :math:`x_{1}(0)=3`, :math:`x_{1}'(0)=0`,
:math:`x_{2}(0)=3`, :math:`x_{2}'(0)=0`.

Solución: Toma la transformada de Laplace de la primera ecuación (con
la notación :math:`x=x_{1}`, :math:`y=x_{2}`):

::

    sage: de1 = maxima("2*diff(x(t),t, 2) + 6*x(t) - 2*y(t)")
    sage: lde1 = de1.laplace("t","s"); lde1.sage()
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

El resultado puede ser difícil de leer, pero significa que

.. math:: -2x'(0) + 2s^2*X(s) - 2sx(0) - 2Y(s) + 6X(s) = 0


(donde la transformada de Laplace de una función en letra minúscula como
:math:`x(t)` es la función en letra mayúscula :math:`X(s)`).
Toma la transformada de Laplace de la segunda ecuación:

::

    sage: t,s = SR.var('t,s')
    sage: x = function('x')
    sage: y = function('y')
    sage: f = 2*x(t).diff(t,2) + 6*x(t) - 2*y(t)
    sage: f.laplace(t,s)
    2*s^2*laplace(x(t), t, s) - 2*s*x(0) + 6*laplace(x(t), t, s) - 2*laplace(y(t), t, s) - 2*D[0](x)(0)

Esto dice

.. math:: -Y'(0) + s^2Y(s) + 2Y(s) - 2X(s) - sy(0) = 0.


Introduce las condiciones iniciales para :math:`x(0)`, :math:`x'(0)`,
:math:`y(0)` y :math:`y'(0)` y resuelve las dos ecuaciones resultantes:

::

    sage: var('s X Y')
    (s, X, Y)
    sage: eqns = [(2*s^2+6)*X-2*Y == 6*s, -2*X +(s^2+2)*Y == 3*s]
    sage: solve(eqns, X,Y)
    [[X == 3*(s^3 + 3*s)/(s^4 + 5*s^2 + 4),
      Y == 3*(s^3 + 5*s)/(s^4 + 5*s^2 + 4)]]

Ahora toma la transformada inversa de Laplace para obtener la respuesta:

::

    sage: var('s t')
    (s, t)
    sage: inverse_laplace((3*s^3 + 9*s)/(s^4 + 5*s^2 + 4),s,t)
    cos(2*t) + 2*cos(t)
    sage: inverse_laplace((3*s^3 + 15*s)/(s^4 + 5*s^2 + 4),s,t)
    -cos(2*t) + 4*cos(t)

Por tanto, la solución es

.. math:: x_1(t) = \cos(2t) + 2\cos(t), \quad x_2(t) = 4\cos(t) - \cos(2t).


La solución puede dibujarse paramétricamente usando

::

    sage: t = var('t')
    sage: P = parametric_plot((cos(2*t) + 2*cos(t), 4*cos(t) - cos(2*t) ),\
    ....: (0, 2*pi), rgbcolor=hue(0.9))
    sage: show(P)

Los componentes individuales pueden dibujarse usando

::

    sage: t = var('t')
    sage: p1 = plot(cos(2*t) + 2*cos(t), 0, 2*pi, rgbcolor=hue(0.3))
    sage: p2 = plot(4*cos(t) - cos(2*t), 0, 2*pi, rgbcolor=hue(0.6))
    sage: show(p1 + p2)

REFERENCIAS: Nagle, Saff, Snider, Fundamentos De Ecuaciones
Diferenciales, 6a ed, Addison-Wesley, 2004. (véase § 5.5).

Método De Euler Para Sistemas De Ecuaciones Diferenciales
---------------------------------------------------------

En el siguiente ejemplo, ilustraremos el método de Euler para EDOs
de primer y segundo órden. Primero, recordemos la idea básica para
ecuaciones de primer órden. Dado un problema con valor inicial de la forma

.. math::

    y'=f(x,y)
    y(a)=c

queremos encontrar el valor aproximado de la solución en :math:`x=b` con :math:`b>a`.

Recuerda de la definición de derivada que

.. math::  y'(x) \approx \frac{y(x+h)-y(x)}{h},


donde :math:`h>0` está dado y es pequeño. Esto, junto con la ED, dan
:math:`f(x,y(x))\approx \frac{y(x+h)-y(x)}{h}`. Ahora resuelve para :math:`y(x+h)`:

.. math::   y(x+h) \approx y(x) + h*f(x,y(x)).


Si llamamos a :math:`h f(x,y(x))` el "término de corrección" (a falta de
algo mejor), llamamos a :math:`y(x)` "el valor viejo de *y*", y
llamamos a :math:`y(x+h)` el "nuevo valor de *y*", entonces, esta
aproximación puede re-expresarse como

.. math::   y_{nuevo} \approx y_{viejo} + h*f(x,y_{viejo}).


Si descomponemos el intervalo desde *a* a *b* en *n* pasos, de modo que
:math:`h=\frac{b-a}{n}`, podemos guardar la información dada por
este método en una tabla.

============== ==================   ================
:math:`x`      :math:`y`            :math:`hf(x,y)`
============== ==================   ================
:math:`a`      :math:`c`            :math:`hf(a,c)`
:math:`a+h`    :math:`c+hf(a,c)`    ...
:math:`a+2h`   ...
...
:math:`b=a+nh` ???                  ...
============== ==================   ================


La meta es llenar todos los espacios de la tabla, una fila cada
la vez, hasta que lleguemos a la casilla ???, que será la
aproximación del método de Euler para :math:`y(b)`.

La idea para los sistemas de EDOs es similar.

**Ejemplo:** Aproxima numéricamente :math:`z(t)` en :math:`t=1` usando 4
pasos del método de Euler, donde :math:`z''+tz'+z=0`,
:math:`z(0)=1`, :math:`z'(0)=0`.

Debemos reducir la EDO de segundo órden a un sistema de dos EDs
de primer órden (usando :math:`x=z`, :math:`y=z'`) y aplicar el método de Euler:

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

Por tanto, :math:`z(1)\approx 0.75`.

También podemos dibujar los puntos :math:`(x,y)` para obtener una representación
aproximada de la curva. La función que hace esto es ``eulers_method_2x2_plot``.
Para poder usarla, necesitamos definir las funciones *f* y
*g* que toman un argumento con tres coordenadas: (*t*, *x*,*y*).

::

    sage: f = lambda z: z[2]        # f(t,x,y) = y
    sage: g = lambda z: -sin(z[1])  # g(t,x,y) = -sin(x)
    sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)

A estas alturas, ``P`` está guardando dos gráficas: ``P[0]``, el gráfico de *x*
vs. *t*, y ``P[1]``, el gráfico de *y* vs. *t*. Podemos mostrar ámbas como sigue:

.. link

::

    sage: show(P[0] + P[1])


Funciones Especiales
--------------------

Se han implementado varios polinomios ortogonales y funciones especiales,
utilizando tanto PARI [GAP]_ como Maxima [Max]_. Estas funciones están
documentadas en las secciones apropiadas ("Polinomios Ortogonales"
y "Funciones Especiales", respectivamente) del manual de referencia de Sage.

::

    sage: x = polygen(QQ, 'x')
    sage: chebyshev_U(2,x)
    4*x^2 - 1
    sage: bessel_I(1,1).n(250)
    0.56515910399248502720769602760986330732889962162109200948029448947925564096
    sage: bessel_I(1,1).n()
    0.565159103992485
    sage: bessel_I(2,1.1).n()  # los últimos digitos son al azar
    0.16708949925104...

Hasta este punto, Sage únicamente ha encapsulado estas funciones para uso numérico.
Para uso simbólico, por favor utiliza directamente la interfaz a Maxima, como en
el siguiente ejemplo:

::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'
    sage: maxima.eval("diff(f,w)")
    '(bessel_y(v-1,w)-bessel_y(v+1,w))/2'

.. [GAP] El Grupo GAP, ``GAP - Grupos, Algorítmos y Programación``, https://www.gap-system.org

.. [Max] Maxima, http://maxima.sf.net/
