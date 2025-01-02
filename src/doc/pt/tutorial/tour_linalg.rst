.. _section-linalg:

Álgebra Linear
==============

O Sage fornece os objetos usuais em álgebra linear, por exemplo, o
polinômio característico, matriz escalonada, traço, decomposição,
etc., de uma matriz.

Criar e multiplicar matrizes é fácil e natural:

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

Note que no Sage, o núcleo de uma matriz :math:`A` é o núcleo à
esquerda, i.e., o conjunto de vetores :math:`w` tal que :math:`wA=0`.

Resolver equações matriciais é fácil usando o método ``solve_right``.
Calculando ``A.solve_right(Y)`` obtém-se uma matrix (ou vetor)
:math:`X` tal que :math:`AX=Y`:

.. link

::

    sage: Y = vector([0, -4, -1])
    sage: X = A.solve_right(Y)
    sage: X
    (-2, 1, 0)
    sage: A * X   # checking our answer...
    (0, -4, -1)

Se não existir solução, o Sage retorna um erro:

.. skip

::

    sage: A.solve_right(w)
    Traceback (most recent call last):
    ...
    ValueError: matrix equation has no solutions

Similarmente, use ``A.solve_left(Y)`` para resolver para :math:`X` em
:math:`XA=Y`.

O Sage também pode calcular autovalores e autovetores::

    sage: A = matrix([[0, 4], [-1, 0]])
    sage: A.eigenvalues ()
    [-2*I, 2*I]
    sage: B = matrix([[1, 3], [3, 1]])
    sage: B.eigenvectors_left()
    [(4, [(1, 1)], 1), (-2, [(1, -1)], 1)]

(A sintaxe para a resposta de ``eigenvectors_left`` é uma lista com
três componentes: (autovalor, autovetor, multiplicidade).) Autovalores
e autovetores sobre ``QQ`` ou ``RR`` também podem ser calculados
usando o Maxima (veja :ref:`section-maxima`).

Como observado em :ref:`section-rings`, o anel sobre o qual a matriz
esta definida afeta alguma de suas propriedades. A seguir, o primeiro
argumento do comando ``matrix`` diz para o Sage considerar a matriz
como uma matriz de inteiros (o caso ``ZZ``), uma matriz de números
racionais (``QQ``), ou uma matriz de números reais (``RR``)::

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

Espaços de Matrizes
-------------------

Agora criamos o espaço :math:`\text{Mat}_{3\times 3}(\QQ)` de matrizes
`3 \times 3` com entradas racionais::

    sage: M = MatrixSpace(QQ,3)
    sage: M
    Full MatrixSpace of 3 by 3 dense matrices over Rational Field

(Para especificar o espaço de matrizes 3 por 4, você usaria
``MatrixSpace(QQ,3,4)``. Se o número de colunas é omitido, ele é
considerado como igual ao número de linhas, portanto,
``MatrixSpace(QQ,3)`` é sinônimo de ``MatrixSpace(QQ,3,3)``.) O espaço
de matrizes é equipado com sua base canônica:

.. link

::

    sage: B = M.basis()
    sage: len(B)
    9
    sage: B[0,1]
    [0 1 0]
    [0 0 0]
    [0 0 0]

Vamos criar uma matriz como um elemento de ``M``.

.. link

::

    sage: A = M(range(9)); A
    [0 1 2]
    [3 4 5]
    [6 7 8]

A seguir calculamos a sua forma escalonada e o núcleo.

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

Agora ilustramos o cálculo com matrizes definidas sobre um corpo
finito:

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

Criamos o subespaço sobre `\GF{2}` gerado pelas linhas acima.

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

A base de `S` usada pelo Sage é obtida a partir das linhas não-nulas
da forma escalonada da matriz de geradores de `S`.

Álgebra Linear Esparsa
----------------------

O Sage fornece suporte para álgebra linear esparsa.

::

    sage: M = MatrixSpace(QQ, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()                  

O algoritmo multi-modular no Sage é bom para matrizes quadradas (mas
não muito bom para matrizes que não são quadradas):

::

    sage: M = MatrixSpace(QQ, 50, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()                  
    sage: M = MatrixSpace(GF(2), 20, 40, sparse=True)
    sage: A = M.random_element()
    sage: E = A.echelon_form()

Note que o Python é sensível a maiúsculas e minúsculas:

::

    sage: M = MatrixSpace(QQ, 10,10, Sparse=True)
    Traceback (most recent call last):
    ...
    TypeError: ...__init__() got an unexpected keyword argument 'Sparse'
