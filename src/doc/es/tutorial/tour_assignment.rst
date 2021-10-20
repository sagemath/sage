
Asignación, Igualdad y Aritmética
====================================

Con algunas excepciones menores, Sage utiliza el lenguaje de programación Python,
de modo que la mayoría de los libros introductorios sobre Python te ayudarán a aprender Sage.

Sage utiliza ``=`` para la asignación. Utiliza ``==``, ``<=``, ``>=``, ``<`` y ``>`` para
la comparación:

::

    sage: a = 5
    sage: a
    5
    sage: 2 == 2
    True
    sage: 2 == 3
    False
    sage: 2 < 3
    True
    sage: a == 5
    True

Sage provee todo lo relacionado con las operaciones matemáticas básicas:

::

    sage: 2**3    #  ** significa exponente
    8
    sage: 2^3     #  ^ es un sinónimo de ** (diferente de Python)
    8
    sage: 10 % 3  #  para argumentos enteros, % significa mod, es decir, resíduo
    1
    sage: 10/4
    5/2
    sage: 10//4   #  para argumentos enteros, // devuelve el cociente de enteros
    2
    sage: 4 * (10 // 4) + 10 % 4 == 10
    True
    sage: 3^2*4 + 2%5
    38

El cálculo de una expresión tal como ``3^2*4 + 2%5`` depende de
el órden en que las operaciones son aplicadas.

Sage también provee muchas funciones matemáticas conocidas; he aquí
solo unos cuantos ejemplos:

::

    sage: sqrt(3.4)
    1.84390889145858
    sage: sin(5.135)
    -0.912021158525540
    sage: sin(pi/3)
    1/2*sqrt(3)

Como demuestra el último ejemplo, algunas expresiones matemáticas devuelven
valores 'exactos', en lugar de aproximaciones numéricas. Para obtener una
aproximación numérica, utilice la función ``n`` o el método
``n`` (ámbas tienen un nombre más largo, ``numerical_approx``, y
la función ``N`` es la misma que ``n``)). Éstas toman argumentos opcionales
``prec``, que es el número requerido de bits de precisión, y ``digits``,
que es el número requerido de digitos decimales de precisión;
el número predeterminado es de 53 bits de precisión.

::

    sage: exp(2)
    e^2
    sage: n(exp(2))
    7.38905609893065
    sage: sqrt(pi).numerical_approx()
    1.77245385090552
    sage: sin(10).n(digits=5)
    -0.54402
    sage: N(sin(10),digits=10)
    -0.5440211109
    sage: numerical_approx(pi, prec=200)
    3.1415926535897932384626433832795028841971693993751058209749

Python es un lenguaje de tipado dinámico, de modo que el valor referido por cada
variable tiene un tipo asociado. Pero una variable dada puede
contener valores de cualquier tipo Python dentro de un ámbito dado:

::

    sage: a = 5   # a es un entero
    sage: type(a)
    <class 'sage.rings.integer.Integer'>
    sage: a = 5/3  # ahora es un número racional
    sage: type(a)
    <class 'sage.rings.rational.Rational'>
    sage: a = 'hello'  # ahora es una cadena
    sage: type(a)
    <... 'str'>

El lenguaje de programación C, que es un lenguaje de tipado estático, es muy
diferente; una variable declarada como int solo puede contener un int
en su ámbito.
