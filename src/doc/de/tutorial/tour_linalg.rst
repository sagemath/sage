.. _section-linalg:

Lineare Algebra
===============

Sage stellt standardmäßige Konstruktionen der Linearen Algebra zur
Verfügung. Zum Beispiel das charakteristische Polynom, die
Zeilenstufenform, die Spur, die Zerlegung von Matrizen, usw..


Das Erzeugen von Matrizen und die Matrixmultiplikation sind einfach
und natürlich:

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

Beachten Sie, dass in Sage der Kern einer Matrix :math:`A` der "linke
Kern", d.h. der Raum der Vektoren :math:`w` mit :math:`wA=0` ist.

Mit der Methode ``solve_right`` können Matrixgleichungen einfach
gelöst werden. Das Auswerten von ``A.solve_right(Y)`` gibt eine Matrix
(oder einen Vektor) :math:`X` zurück, so dass :math:`AX=Y` gilt:

.. link

::

    sage: Y = vector([0, -4, -1])
    sage: X = A.solve_right(Y)
    sage: X
    (-2, 1, 0)
    sage: A * X   # wir überprüfen unsere Antwort...
    (0, -4, -1)

Falls keine Lösung existiert, gibt Sage einen Fehler zurück:

.. skip

::

    sage: A.solve_right(w)
    Traceback (most recent call last):
    ...
    ValueError: matrix equation has no solutions

Auf ähnliche Weisen können Sie ``A.solve_left(Y)`` benutzen um nach :math:`X` in
:math:`XA=Y` aufzulösen.

Sage kann auch Eigenwerte und Eigenvektoren berechnen::

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

(Die Syntax der Ausgabe von ``eigenvectors_left`` ist eine Liste von
Tripeln: (Eigenwert, Eigenvektor, Vielfachheit).) Eigenwerte und
Eigenvektoren über ``QQ`` oder ``RR`` können auch unter Verwendung von
Maxima berechnen werden (Lesen Sie :ref:`section-maxima` unterhalb).

Wie in :ref:`section-rings` bemerkt wurde, beeinflusst der Ring, über
dem die Matrix definiert ist, einige ihrer Eigenschaften. Im Folgenden
gibt erste Argument des ``matrix``-Befehls Sage zu verstehen, dass die
Matrix als Matrix über den ganzen Zahlen (``ZZ``), als Matrix über den
rationalen Zahlen (``QQ``), oder als Matrix über den reellen Zahlen
(``RR``), aufgefasst werden soll::

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

Um Eigenwerte und Eigenvektoren mit reellen oder komplexen Gleitkommazahlen zu
berechnen sollte die Matrix über ``RDF`` (Real Double Field = Körper der
reellen Gleitkommazahlen mit doppelter Genauigkeit) oder ``CDF`` (Complex Double
Field = Körper der komplexen Gleitkommazahlen mit doppelter Genauigkeit)
definiert werden. Falls kein Koeffizientenring angegeben wird und die
Matrixeinträge relle oder komplexe Gleitkommazahlen sind dann werden
standardmässig die Körper ``RR`` oder ``CC`` verwendet, welche allerdings nicht
alle der folgenden Berechnungen unterstützen::

    sage: ARDF = matrix(RDF, [[1.2, 2], [2, 3]])
    sage: ARDF.eigenvalues()  # rel tol 8e-16
    [-0.09317121994613098, 4.293171219946131]
    sage: ACDF = matrix(CDF, [[1.2, I], [2, 3]])
    sage: ACDF.eigenvectors_right()  # rel tol 3e-15
    [(0.8818456983293743 - 0.8209140653434135*I, [(0.7505608183809549, -0.616145932704589 + 0.2387941530333261*I)], 1),
    (3.3181543016706256 + 0.8209140653434133*I, [(0.14559469829270957 + 0.3756690858502104*I, 0.9152458258662108)], 1)]


Matrizenräume
-------------

Wir erzeugen den Raum :math:`\text{Mat}_{3\times 3}(\QQ)` der  `3 \times
3` Matrizen mit rationalen Einträgen::

    sage: M = MatrixSpace(QQ,3)
    sage: M
    Full MatrixSpace of 3 by 3 dense matrices over Rational Field

(Um den Raum der 3 mal 4 Matrizen anzugeben würden Sie
``MatrixSpace(QQ,3,4)`` benutzen. Falls die Anzahl der Spalten nicht
angegeben wurde, ist diese standardmäßig gleich der Anzahl der Zeilen,
so dass ``MatrixSpace(QQ,3)`` ein Synonym für ``MatrixSpace(QQ,3,3)``
ist.) Der Matrizenraum ist mit seiner kanonischen Basis ausgestattet:

.. link

::

    sage: B = M.basis()
    sage: len(B)
    9
    sage: B[0,1]
    [0 1 0]
    [0 0 0]
    [0 0 0]

Wir erzeugen eine Matrix als ein Element von ``M``.

.. link

::

    sage: A = M(range(9)); A
    [0 1 2]
    [3 4 5]
    [6 7 8]

Als nächstes berechnen wir die reduzierte Zeilenstufenform und den Kern.

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

Nun zeigen wir, wie man Matrizen berechnen, die über
endlichen Körpern definiert sind:

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

Wir erstellen den Unterraum von `\GF{2}^8`, der von den obigen Zeilen
aufgespannt wird.

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

Die Basis von `S`, die von Sage benutzt wird, wird aus den von Null
verschiedenen Zeilen der reduzierten Zeilenstufenform der Matrix der
Generatoren von `S` erhalten.

Lineare Algebra mit dünnbesetzten Matrizen
------------------------------------------

Sage unterstützt Lineare Algebra mit dünnbesetzten Matrizen über
Hauptidealringen.


::

    sage: M = MatrixSpace(QQ, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()

Der multi-modulare Algorithmus kann bei quadratischen Matrizen gut
angewendet werden (bei nicht quadratischen Matrizen ist er nicht so gut):

::

    sage: M = MatrixSpace(QQ, 50, 100, sparse=True)
    sage: A = M.random_element(density = 0.05)
    sage: E = A.echelon_form()
    sage: M = MatrixSpace(GF(2), 20, 40, sparse=True)
    sage: A = M.random_element()
    sage: E = A.echelon_form()

Beachten Sie, dass Python zwischen Klein- und Großschreibung unterscheidet:

::

    sage: M = MatrixSpace(QQ, 10,10, Sparse=True)
    Traceback (most recent call last):
    ...
    TypeError: ...__init__() got an unexpected keyword argument 'Sparse'...
