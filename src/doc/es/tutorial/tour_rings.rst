.. -*- coding: utf-8 -*-
.. _section-rings:

Anillos Elementales
===================

Cuando definimos matrices, vectores o polinomios, a veces es útil,
incluso necesario, especificar el "anillo" sobre el que están definidos.
Un *anillo* es una construcción matemática consistente en un conjunto de
elementos sobre los que está bien definidas las operaciones de suma y producto;
si la noción de anillo no te resulta familiar, probablemente sólo necesitas 
conocer estos cuatro anillos:

* los enteros `\{..., -1, 0, 1, 2, ...\}`, a los que nos referimos en Sage 
  por ``ZZ``.
* los números racionales -- e.g., fracciones, o cocientes de números enteros 
  --, ``QQ`` en Sage.
* los números reales, ``RR`` en Sage.
* los números complejos, ``CC`` en Sage.

Es importante conocer estas distinciones porque el mismo polinomio, por 
ejemplo, puede ser tratado de forma diferente dependiendo del anillo sobre el
que se ha definido. Por ejemplo, el polinomio `x^2-2` tiene dos raíces,
`\pm \sqrt{2}`.  Estas raíces no son racionales, así que si trabajamos
con polinomios con coeficientes racionales, el polinomio es irreducible.
Sin embargo, si los coeficientes son números reales, el polinomio factoriza 
como producto de dos factores lineales. En el siguiente ejemplo, los conjuntos
de polinomios se llaman "ratpoly" y "realpoly", aunque no usaremos estos 
nombres; observa sin embargo que las cadenas ".<t>" y ".<z>" sirven para dar
nombre a las variables usadas en cada caso. ::

    sage: ratpoly.<t> = PolynomialRing(QQ)
    sage: realpoly.<z> = PolynomialRing(RR)

Veamos el punto que hicimos antes sobre factorizar `x^2-2`:

.. link

::

    sage: factor(t^2-2)
    t^2 - 2
    sage: factor(z^2-2)
    (z - 1.41421356237310) * (z + 1.41421356237310)

Comentarios similares se aplican a las matrices: la forma reducida por filas
de una matriz puede depender del anillo en que está definida, al igual que
sus autovalores y autofunciones. Hay más construcciones con polinomios en la 
sección :ref:`section-poly`, y más construcciones con matrices en
:ref:`section-linalg`.

El símbolo ``I`` representa la raíz cuadrada de :math:`-1`; ``i`` es un
sinónimo de ``I``. Por supuesto, no es un número racional::

    sage: i  # raíz cuadrada de -1
    I     
    sage: i in QQ
    False

Nota: El código siguiente puede no funcionar como esperas si hemos asignado
otro valor a la variable ``i``, por ejemplo si la hemos usado como variable
interna de un bucle. En este caso, podemos devolver ``i`` a su valor original::

    sage: reset('i')


Hay una sutileza al definir números complejos: el símbolo ``i`` representa una
raíz cuadrada de `-1`, pero es una raíz *formal* como elemento de un cuerpo
de números algebráicos.
Ejecutando ``CC(i)`` o ``CC.0``  o ``CC.gen(0)`` obtenemos el número *complejo*
de coma flotante que es la raíz cuadrada de `-1`. ::

    sage: i = CC(i)       # número complejo de coma flotante
    sage: i == CC.0
    True
    sage: a, b = 4/3, 2/3
    sage: z = a + b*i
    sage: z
    1.33333333333333 + 0.666666666666667*I
    sage: z.imag()        # parte imaginaria
    0.666666666666667
    sage: z.real() == a   # conversión automática antes de la comparación
    True
    sage: a + b
    2
    sage: 2*b == a
    True
    sage: parent(2/3)
    Rational Field
    sage: parent(4/2)
    Rational Field
    sage: 2/3 + 0.1       # conversión automática antes de la suma
    0.766666666666667
    sage: 0.1 + 2/3       # las reglas de conversión son simétricas en Sage
    0.766666666666667

Veamos más ejemplos de anillos elementales en Sage. Como mencionamos antes,
nos podemos referir al anillo de números racionales usando ``QQ``, o también 
``RationalField()`` (*field*, o *cuerpo*, se refiere a un anillo en el que
el producto es conmutativo y todo elemento excepto el cero tiene un inverso 
para la multiplicación. De este modo, los racionales son un cuerpo, pero los
enteros no::

    sage: RationalField()
    Rational Field
    sage: QQ
    Rational Field
    sage: 1/2 in QQ
    True

El número decimal ``1.2`` se considera que está en ``QQ``: los números 
decimales, que también son racionales, se pueden convertir a racionales de
forma automática. Sin embargo, los números `\pi` y `\sqrt{2}` no son 
racionales::

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

En Sage también podemos trabajar con otros anillos, como cuerpos finitos,
enteros `p`-ádicos, el anillo de los números algebraicos, anillos de polinomios 
y anillos de matrices. Veamos algunos de estos anillos::

    sage: GF(3)
    Finite Field of size 3
    sage:                 # es necesario dar un nombre al generador si el número
    sage: GF(27, 'a')     # de elementos no es primo
    Finite Field in a of size 3^3
    sage: Zp(5)
    5-adic Ring with capped relative precision 20
    sage: sqrt(3) in QQbar # clausura algebraica de QQ
    True
