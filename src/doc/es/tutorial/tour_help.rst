.. _chapter-help:

Obteniendo Ayuda
================

Sage posee una extensa documentación incorporada, accesible con solo teclear
el nombre de una función o una constante (por ejemplo), seguido de un signo
de interrogación:

.. skip

::

    sage: tan?
    Type:        <class 'sage.calculus.calculus.Function_tan'>
    Definition:  tan( [noargspec] )
    Docstring:

        The tangent function

        EXAMPLES:
            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
    sage: log2?
    Type:        <class 'sage.functions.constants.Log2'>
    Definition:  log2( [noargspec] )
    Docstring:

        The natural logarithm of the real number 2.

        EXAMPLES:
            sage: log2
            log2
            sage: float(log2)
            0.69314718055994529
            sage: RR(log2)
            0.693147180559945
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(log2)
            0.69314718055994530941723212145817656807550013436025525412068
            sage: l = (1-log2)/(1+log2); l
            (1 - log(2))/(log(2) + 1)
            sage: R(l)
            0.18123221829928249948761381864650311423330609774776013488056
            sage: maxima(log2)
            log(2)
            sage: maxima(log2).float()
            .6931471805599453
            sage: gp(log2)
            0.6931471805599453094172321215             # 32-bit
            0.69314718055994530941723212145817656807   # 64-bit
    sage: sudoku?
    File:        sage/local/lib/python2.5/site-packages/sage/games/sudoku.py
    Type:        <... 'function'>
    Definition:  sudoku(A)
    Docstring:

        Solve the 9x9 Sudoku puzzle defined by the matrix A.

        EXAMPLE:
            sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0,
        0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0, 0,0,0, 2,0,8, 0,0,0,
        0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2,
        0,0,0, 4,9,0, 0,5,0, 0,0,3])
            sage: A
            [5 0 0 0 8 0 0 4 9]
            [0 0 0 5 0 0 0 3 0]
            [0 6 7 3 0 0 0 0 1]
            [1 5 0 0 0 0 0 0 0]
            [0 0 0 2 0 8 0 0 0]
            [0 0 0 0 0 0 0 1 8]
            [7 0 0 0 0 4 1 5 0]
            [0 3 0 0 0 2 0 0 0]
            [4 9 0 0 5 0 0 0 3]
            sage: sudoku(A)
            [5 1 3 6 8 7 2 4 9]
            [8 4 9 5 2 1 6 3 7]
            [2 6 7 3 4 9 5 8 1]
            [1 5 8 4 6 3 9 7 2]
            [9 7 4 2 1 8 3 6 5]
            [3 2 6 7 9 5 4 1 8]
            [7 8 2 9 3 4 1 5 6]
            [6 3 5 1 7 2 8 9 4]
            [4 9 1 8 5 6 7 2 3]

Sage también provee 'Autocompletado con el tabulador': teclea las primeras letras de
una función y luego oprime la tecla del tabulador. Por ejemplo, si tecleas ``ta``
seguido por ``TAB``, Sage imprimirá
``tachyon, tan, tanh, taylor``.
Esto proporciona una buena manera de encontrar los nombres de funciones y otras
estructuras en Sage.


.. _section-functions:

Funciones, Indentación Y Conteo
====================================

Para definir una nueva función en Sage, utilice el comando ``def`` y el signo de dos puntos
después de la lista de nombres de variable. Por ejemplo:

::

    sage: def is_even(n):
    ....:     return n % 2 == 0
    ....:
    sage: is_even(2)
    True
    sage: is_even(3)
    False

Nota: Dependiendo de la versión del tutorial que estás leyendo, puede que veas 
puntos ``....:`` en la segunda línea de este ejemplo.
No los incluyas; son solo para enfatizar que el código está indentado.
Siempre que este sea el caso, presiona [Return/Enter] una vez al final del bloque
para insertar una línea en blanco y concluir la definición de la función.

No tienes que especificar los tipos de ninguno de los argumentos de entrada.
Puedes especificar múltiples entradas, cada una de las cuales puede tener un
valor predeterminado opcional. Por ejemplo, la función de abajo tiene un valor
predeterminado ``divisor=2`` si no se especifica el valor de ``divisor``.

::

    sage: def is_divisible_by(number, divisor=2):
    ....:     return number % divisor == 0
    sage: is_divisible_by(6,2)
    True
    sage: is_divisible_by(6)
    True
    sage: is_divisible_by(6, 5)
    False

También puedes especificar explícitamente una o ambas de las entradas cuando
llames a la función; si especificas las entradas explícitamente, puedes darlas
en cualquier órden:

.. link

::

    sage: is_divisible_by(6, divisor=5)
    False
    sage: is_divisible_by(divisor=2, number=6)
    True

En Python, los bloques de código no se encierran entre llaves o bloques begin...end
como en muchos otros lenguajes. En vez de ello, los bloques de código
se indican por medio de la indentación, la cual se debe agrupar con exactitud.
Por ejemplo, el siguiente es un error de sintáxis porque la declaración ``return``
no está indentada al mismo nivel que las otras líneas por encima de ella.

.. skip

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:    return v
    Syntax Error:
           return v

Si arreglas la indentación, la función se ejecutará:

::

    sage: def even(n):
    ....:     v = []
    ....:     for i in range(3,n):
    ....:         if i % 2 == 0:
    ....:             v.append(i)
    ....:     return v
    sage: even(10)
    [4, 6, 8]

El punto y coma no es necesario al final de las líneas. Una línea termina, en muchos casos,
por un carácter de nueva línea. Sin embargo, puedes poner múltiples declaraciones
en una línea, separadas por punto y coma:

::

    sage: a = 5; b = a + 3; c = b^2; c
    64

Si quisieras que una simple línea de código abarque multiples líneas, utiliza
una barra invertida como terminación:

::

    sage: 2 + \
    ....: 3
    5

En Sage, se cuenta iterando sobre un rango de enteros. Por ejemplo,
la primer línea de abajo es exactamente igual a ``for(i=0; i<3; i++)`` en C++ o Java:

::

    sage: for i in range(3):
    ....:     print(i)
    0
    1
    2

La primer línea de abajo es igual a ``for(i=2;i<5;i++)``.

::

    sage: for i in range(2,5):
    ....:     print(i)
    2
    3
    4

El tercer argumento controla el incremento, de modo que lo siguiente es igual a
``for(i=1;i<6;i+=2)``.

::

    sage: for i in range(1,6,2):
    ....:     print(i)
    1
    3
    5

A menudo, querrás crear una tabla para presentar números que has calculado
utilizando Sage. Una manera sencilla de hacer esto es usando el formateado de
cadenas. Abajo, creamos tres columnas, cada una con un ancho exácto de 6 caracteres 
y hacemos una tabla de cuadrados y cubos.

::

    sage: for i in range(5):
    ....:     print('%6s %6s %6s' % (i, i^2, i^3))
         0      0      0
         1      1      1
         2      4      8
         3      9     27
         4     16     64

La estructura de datos más básica en Sage es la lista, la cual es -- como
sugiere su nombre -- solo una lista de objetos arbitrarios.
Por ejemplo, el comando ``range`` que hemos usado crea una lista::

    sage: list(range(2,10))
    [2, 3, 4, 5, 6, 7, 8, 9]

He aquí una lista más complicada:

::

    sage: v = [1, "hello", 2/3, sin(x^3)]
    sage: v
    [1, 'hello', 2/3, sin(x^3)]

El indexado de una lista comienza en el cero, como en muchos lenguajes de programación.

.. link

::

    sage: v[0]
    1
    sage: v[3]
    sin(x^3)

La función ``len(v)`` devuelve la longitud de ``v``. Utiliza ``v.append(obj)`` para
añadir un nuevo objeto al final de ``v``, y utiliza ``del v[i]`` para borrar
el :math:`i-ésimo` elemento de ``v``:

.. link

::

    sage: len(v)
    4
    sage: v.append(1.5)
    sage: v
    [1, 'hello', 2/3, sin(x^3), 1.50000000000000]
    sage: del v[1]
    sage: v
    [1, 2/3, sin(x^3), 1.50000000000000]

Otra estructura de datos importante es el diccionario (o array asociativo).
Funciona como una lista, excepto que puede ser indexado con casi
cualquier objeto (los índices deben ser immutables):

::

    sage: d = {'hi':-2,  3/8:pi,   e:pi}
    sage: d['hi']
    -2
    sage: d[e]
    pi

También puedes definir nuevos tipos de datos usando clases. El encapsulado
de objetos matemáticos con clases es una técnica potente que puede
ayudar a simplificar y organizar tus programas en Sage. Abajo, definimos una
clase que representa la lista de enteros positivos pares hasta ``n``;
se deriva de el tipo básico ``list``.

::

    sage: class Evens(list):
    ....:     def __init__(self, n):
    ....:         self.n = n
    ....:         list.__init__(self, range(2, n+1, 2))
    ....:     def __repr__(self):
    ....:         return "Even positive numbers up to n."

El método ``__init__`` se llama para inicializar al objeto cuando
es creado; el método ``__repr__`` imprime el objeto.
Llamamos al método constructor de listas en la segunda línea del
método ``__init__``. A continuación, creamos un objeto de clase ``Evens``:

.. link

::

    sage: e = Evens(10)
    sage: e
    Even positive numbers up to n.

Observe que ``e`` se imprime usando el método ``__repr__`` que hemos definido.
Para ver la lista subyacente de números, utilice la función ``list``:

.. link

::

    sage: list(e)
    [2, 4, 6, 8, 10]

También podemos acceder al atributo ``n`` o tratar a ``e`` como una lista.

.. link

::

    sage: e.n
    10
    sage: e[2]
    6
