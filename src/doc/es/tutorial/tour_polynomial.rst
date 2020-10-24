.. -*- coding: utf-8 -*-
.. _section-poly:

Polinomios
===========

En esta sección mostraremos cómo crear y usar polinomios en Sage.


.. _section-univariate:

Polinomios en una variable
--------------------------

Hay tres formas de construir anillos de polinomios.

::

    sage: R = PolynomialRing(QQ, 't')
    sage: R
    Univariate Polynomial Ring in t over Rational Field

De esta forma creamos un anillo de polinomios en una variable, y pedimos
que esta variable se muestre por pantalla como ``t``. Sin embargo, de esta forma
no se define ``t`` como variable simbólica en Sage, y no se puede usar este
símbolo para escribr polinomios de ``R``  como por ejemplo  :math:`t^2+1`.

Otra forma es:

.. link

::

    sage: S = QQ['t']
    sage: S == R
    True

Los mismos comentarios sobre la variable ``t`` se aplican a esta forma.

Una tercera forma, muy práctica es

::

    sage: R.<t> = PolynomialRing(QQ)

o

::

    sage: R.<t> = QQ['t']

o incluso

::

    sage: R.<t> = QQ[]

Todas estas formas tienen el efecto añadido de definir la variable ``t`` como
la indeterminada del anillo de polinomios, lo que hace más sencillo definir
elementos de ``R``. (Esta tercera forma es similar a la notación
de Magma [MAGMA]_ , y al igual que en Magma se puede usar para una amplia variedad de
objetos.)

.. link

::

    sage: poly = (t+1) * (t+2); poly
    t^2 + 3*t + 2
    sage: poly in R
    True

Independientemente de la forma usada para definir un anillo de polinomios, 
podemos recuperar la indeterminada mediante el generador :math:`0`-ésimo.

::

    sage: R = PolynomialRing(QQ, 't')
    sage: t = R.0
    sage: t in R
    True

Observa que una construcción similar funciona con los números complejos, que
pueden ser vistos como el conjunto generado por los números reales y el
símbolo ``i``:

::

    sage: CC
    Complex Field with 53 bits of precision
    sage: CC.0  # 0th generator of CC
    1.00000000000000*I

También podemos obtener tanto el anillo como el generador, o sólo el generador, 
en el momento de crear un anillo de polinomios, del modo siguiente:

::

    sage: R, t = QQ['t'].objgen()
    sage: t    = QQ['t'].gen()
    sage: R, t = objgen(QQ['t'])
    sage: t    = gen(QQ['t'])

Finalmente hacemos un poco de aritmética en :math:`\QQ[t]`.

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

Observamos que la factorización tiene en cuenta la unidad que multiplica a los
factores irreducibles.

Si en el curso de nuestra investigación usásemos mucho, por ejemplo, la función
``R.cyclotomic_polynomial``, sería recomendable citar, además de a Sage,
a la componente de Sage que realiza el cálculo en última instancia.
En este caso, ejecutando ``R.cyclotomic_polynomial??`` para ver el código 
fuente, observamos la línea ``f = pari.polcyclo(n)`` , lo que significa que
para este cálculo se usa PARI, y deberíamos citarlo además de Sage.

Al dividir dos polinomios, construimos un elemento del cuerpo de fracciones 
(que Sage crea automáticamente).

::

    sage: x = QQ['x'].0
    sage: f = x^3 + 1; g = x^2 - 17
    sage: h = f/g;  h
    (x^3 + 1)/(x^2 - 17)
    sage: h.parent()
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

Usando series de Laurent, podemos calcular expansiones en serie de potencias
de elementos del cuerpo de fracciones de ``QQ[x]``:

::

    sage: R.<x> = LaurentSeriesRing(QQ); R
    Laurent Series Ring in x over Rational Field
    sage: 1/(1-x) + O(x^10)
    1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)

Si usamos otro nombre para la variable, obtenemos un anillo diferente.

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

El anillo de polinomios está determinado por el anillo de coeficientes y la
variable. Observamos que construir otro anillo con una variable de nombre 
``x`` no devuelve un anillo distinto.

::

    sage: R = PolynomialRing(QQ, "x")
    sage: T = PolynomialRing(QQ, "x")
    sage: R == T
    True      
    sage: R is T
    True
    sage: R.0 == T.0
    True

Sage soporta los anillos de series de potencias y de series de Laurent sobre
cualquier anillo base. En el ejemplo siguiente, creamos un elemento de
:math:`\GF{7}[[T]]` y calculamos su inverso para crear un elemento de
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

También podemos crear anillos de series de potencias usando dobles corchetes:

::

    sage: GF(7)[['T']]
    Power Series Ring in T over Finite Field of size 7

Polinomios en varias variables
------------------------------

Para trabajar con polinomios de varias variables, comenzamos por declarar el 
anillo de polinomios y las variables.

::

    sage: R = PolynomialRing(GF(5),3,"z") # here, 3 = number of variables
    sage: R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Al igual que al definir anillos de polinomios en una variable, hay varias 
formas:

::

    sage: GF(5)['z0, z1, z2']
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5
    sage: R.<z0,z1,z2> = GF(5)[]; R
    Multivariate Polynomial Ring in z0, z1, z2 over Finite Field of size 5

Es posible usar una letra distinta para cada variable usando la notación:

::

    sage: PolynomialRing(GF(5), 3, 'xyz')
    Multivariate Polynomial Ring in x, y, z over Finite Field of size 5

Veamos un poco de aritmética:

::

    sage: z = GF(5)['z0, z1, z2'].gens()
    sage: z
    (z0, z1, z2)
    sage: (z[0]+z[1]+z[2])^2
    z0^2 + 2*z0*z1 + z1^2 + 2*z0*z2 + 2*z1*z2 + z2^2

Es posible usar una notación más parecida a la convención usual en matemáticas
para definir el anillo.

::

    sage: R = GF(5)['x,y,z']
    sage: x,y,z = R.gens()
    sage: QQ['x']
    Univariate Polynomial Ring in x over Rational Field
    sage: QQ['x,y'].gens()
    (x, y)
    sage: QQ['x'].objgens()
    (Univariate Polynomial Ring in x over Rational Field, (x,))

Los polinomios en varias variables están implementados en Sage usando 
diccionarios de Python y la "representación distributiva" de un polinomio.
Sage usa en parte Singular [Si]_, por ejemplo para el cálculo del mcd de dos
polinomios y la base de Gröbner de un ideal.

::

    sage: R, (x, y) = PolynomialRing(RationalField(), 2, 'xy').objgens()
    sage: f = (x^3 + 2*y^2*x)^2
    sage: g = x^2*y^2
    sage: f.gcd(g)
    x^2

A continuación creamos el ideal :math:`(f,g)` generado por :math:`f` y
:math:`g`, simplemente multiplicando la tupla ``(f,g)`` por ``R`` (también
podemos escribir ``ideal([f,g])`` o ``ideal(f,g)``).

.. link

::

    sage: I = (f, g)*R; I
    Ideal (x^6 + 4*x^4*y^2 + 4*x^2*y^4, x^2*y^2) of Multivariate Polynomial 
    Ring in x, y over Rational Field
    sage: B = I.groebner_basis(); B
    [x^6, x^2*y^2]
    sage: x^2 in I
    False

La base de Gröbner de arriba no es una lista, sino una secuencia inmutable. 
Esto implica que tiene un universo y un padre, y que no se puede cambiar
(lo cual es importante porque otras rutinas usarán esta base de Gröbner).

.. link

::

    sage: B.parent()
    <class 'sage.rings.polynomial.multi_polynomial_sequence.PolynomialSequence_generic'>
    sage: B.universe()
    Multivariate Polynomial Ring in x, y over Rational Field
    sage: B[1] = x
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.

Sage incluye código basado en la librería Singular que permite hacer algo de 
álgebra conmutativa  (entiéndase: no tanta como nos gustaría). Por ejemplo, 
podemos calcular la descomposición primaria y los primos asociados a :math:`I`:

.. link

::

    sage: I.primary_decomposition()
    [Ideal (x^2) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y^2, x^6) of Multivariate Polynomial Ring in x, y over Rational Field]
    sage: I.associated_primes()
    [Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field,
     Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field]

.. [Si] Singular es un sistema de álgebra computerizado para cálculos con
        polinomios, http://www.singular.uni-kl.de

.. [MAGMA] Sistema de algebra computacional, http://magma.maths.usyd.edu.au/magma/
