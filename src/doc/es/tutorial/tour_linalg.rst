.. -*- coding: utf-8 -*-
.. _section-linalg:

Álgebra Lineal
==============

Sage soporta construcciones estándar de álgebra lineal, como el
polinomio característico, la forma escalonada, la traza, 
descomposición, etcétera de una matriz.

La creación de matrices y la multiplicación es sencilla y natural:

::

    sage: A = Matrix([[1,2,3],[3,2,1],[1,1,1]])
    sage: w = vector([1,1,-4])
    sage: w*A
    (0, 0, 0)
    sage: A*w
    (-9, 1, -2)
    sage: kernel(A)
    Free module of degree 3 and rank 1 over Integer Ring
    Echelon basis matrix:
    [ 1  1 -4]

La descripción de ``kernel(A)`` indica que se trata de un
subespacio de dimensión 1 ("rank 1") de un espacio de dimensión 3
("degree 3"). Por el momento, tanto ``kernel(A)`` como el espacio 
ambiente admiten coeficientes enteros ("over Integer Ring").
Finalmente, Sage nos muestra una base escalonada ("Echelon basis").

Observa que en Sage, el núcleo de la matriz :math:`A` es el "núcleo por
la izquierda", e.g. el subespacio formado por los vectores :math:`w` 
tales que :math:`wA=0`.

Resolver ecuaciones matriciales es sencillo, usando el método 
``solve_right`` (resolver por la derecha). Al evaluar 
``A.solve_right(Y)`` obtenemos una matriz (o un vector)
:math:`X` tal que :math:`AX=Y`:

.. link

::

    sage: Y = vector([0, -4, -1])
    sage: X = A.solve_right(Y)
    sage: X
    (-2, 1, 0)
    sage: A * X   # comprobando la solución...
    (0, -4, -1)

Si no hay solución, Sage lanza un error:

.. skip

::

    sage: A.solve_right(w)
    Traceback (most recent call last):
    ...
    ValueError: matrix equation has no solutions

De forma similar, usamos ``A.solve_left(Y)`` para despejar :math:`X` de
la ecuación :math:`XA=Y`.

Sage también puede calcular autovalores ("eigenvalues") y autovectores
("eigenvectors")::

    sage: A = matrix([[0, 4], [-1, 0]])
    sage: A.eigenvalues ()
    [-2*I, 2*I]
    sage: B = matrix([[1, 3], [3, 1]])
    sage: B.eigenvectors_left()
    [(4, [
    (1, 1)
    ], 1), (-2, [
    (1, -1)
    ], 1)]

(La sintaxis de la salida de ``eigenvectors_left`` es una lista de
tuplas: (autovalor, autovector, multiplicidad).)  Los autovalores
y autovectores sobre ``QQ`` o ``RR`` también se pueden calcular
usando Maxima.

Como ya indicamos en :ref:`section-rings`, el anillo sobre el que se 
define una matriz afecta algunas de sus propiedades. En las líneas que 
siguen, el primer argumento al comando ``matrix`` le dice a Sage que
considere la matriz como una matriz de enteros (si el argumento es
``ZZ``), de números racionales (si es ``QQ``), o de números reales 
(si es ``RR``)::

    sage: AZ = matrix(ZZ, [[2,0], [0,1]])
    sage: AQ = matrix(QQ, [[2,0], [0,1]])
    sage: AR = matrix(RR, [[2,0], [0,1]])
    sage: AZ.echelon_form()
    [2 0]
    [0 1]
    sage: AQ.echelon_form()
    [1 0]
    [0 1]
    sage: AR.echelon_form()
    [ 1.00000000000000 0.000000000000000]
    [0.000000000000000  1.00000000000000]

(El comando ``echelon_form`` devuelve una forma escalonada de la matriz)

Espacios de matrices
--------------------

Creamos el espacio :math:`\text{Mat}_{3\times 3}(\QQ)` matrices 
`3 \times 3` con coeficientes racionales::

    sage: M = MatrixSpace(QQ,3)
    sage: M
    Full MatrixSpace of 3 by 3 dense matrices over Rational Field

(Para especificar el espacio de matrices 3 por 4, usaríamos
``MatrixSpace(QQ,3,4)``. Si se omite el número de columnas, se adopta
por defecto el número de filas, de modo que ``MatrixSpace(QQ,3)``
es un sinónimo de ``MatrixSpace(QQ,3,3)``.) El espacio de matrices
está equipado con su base canónica:

.. link

::

    sage: B = M.basis()
    sage: len(B)
    9
    sage: B[0,1]
    [0 1 0]
    [0 0 0]
    [0 0 0]

Creamos una matriz como un elemento de ``M``.

.. link

::

    sage: A = M(range(9)); A
    [0 1 2]
    [3 4 5]
    [6 7 8]

Calculamos su forma escalonada por filas y su núcleo.

.. link

::

    sage: A.echelon_form()
    [ 1  0 -1]
    [ 0  1  2]
    [ 0  0  0]
    sage: A.kernel()
    Vector space of degree 3 and dimension 1 over Rational Field
    Basis matrix:
    [ 1 -2  1]

Ilustramos un cálculo de matrices definidas sobre cuerpos finitos:

::

    sage: M = MatrixSpace(GF(2),4,8)
    sage: A = M([1,1,0,0, 1,1,1,1, 0,1,0,0, 1,0,1,1, 
    ....:        0,0,1,0, 1,1,0,1, 0,0,1,1, 1,1,1,0])
    sage: A
    [1 1 0 0 1 1 1 1]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 1 1 1 1 1 0]
    sage: rows = A.rows()
    sage: A.columns()
    [(1, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 1), (0, 0, 0, 1), 
     (1, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0)]
    sage: rows
    [(1, 1, 0, 0, 1, 1, 1, 1), (0, 1, 0, 0, 1, 0, 1, 1), 
     (0, 0, 1, 0, 1, 1, 0, 1), (0, 0, 1, 1, 1, 1, 1, 0)]

Construimos el subespacio sobre `\GF{2}` engendrado por las filas de 
arriba.

.. link

::

    sage: V = VectorSpace(GF(2),8)
    sage: S = V.subspace(rows)
    sage: S
    Vector space of degree 8 and dimension 4 over Finite Field of size 2
    Basis matrix:
    [1 0 0 0 0 1 0 0]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 0 1 0 0 1 1]
    sage: A.echelon_form()
    [1 0 0 0 0 1 0 0]
    [0 1 0 0 1 0 1 1]
    [0 0 1 0 1 1 0 1]
    [0 0 0 1 0 0 1 1]

La base de `S` usada por Sage se obtiene de las filas no nulas de la
forma escalonada reducida de la matriz compuesta por los generadores
de `S`.

Álgebra Lineal Dispersa
-----------------------

Sage soporta espacios de matrices sobre DIPs almacenados de forma
dispersa.

::

    sage: M = MatrixSpace(QQ, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()                  

El algoritmo multi-modular de Sage es bueno para matrices cuadradas
(pero no tan bueno para matrices no cuadradas):

::

    sage: M = MatrixSpace(QQ, 50, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()                  
    sage: M = MatrixSpace(GF(2), 20, 40, sparse=True)
    sage: A = M.random_element()
    sage: E = A.echelon_form()

