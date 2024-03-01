************
Introducción
************

Completar este tutorial debería llevarte unas 3 o 4 horas. Puedes leerlo en versión HTML o PDF, o desde el
notebook (interfaz interactiva vía web) de Sage (haz click en ``Help``, luego haz click en ``Tutorial`` para trabajar interactivamente en el tutorial desde dentro de Sage).

Aunque gran parte de Sage está implementada usando el lenguaje de programación
Python, no es necesario ningún conocimiento previo de Python para poder leer este tutorial.
En algún punto seguramente querrás aprender Python (¡un lenguaje muy divertido!), y hay muchos
recursos gratuitos excelentes para hacerlo, La Guía Para Principiantes De Python [PyB]_
enumera muchas opciones.
Si tan solo quieres experimentar ligeramente con Sage, este tutorial es el
lugar justo para empezar. Por ejemplo:

::

    sage: 2 + 2
    4
    sage: factor(-2007)
    -1 * 3^2 * 223

    sage: A = matrix(4,4, range(16)); A
    [ 0  1  2  3]
    [ 4  5  6  7]
    [ 8  9 10 11]
    [12 13 14 15]

    sage: factor(A.charpoly())
    x^2 * (x^2 - 30*x - 80)

    sage: m = matrix(ZZ,2, range(4))
    sage: m[0,0] = m[0,0] - 3
    sage: m
    [-3  1]
    [ 2  3]

    sage: E = EllipticCurve([1,2,3,4,5]);
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
    over Rational Field
    sage: E.anlist(10)
    [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    sage: E.rank()
    1

    sage: k = 1/(sqrt(3)*I + 3/4 + sqrt(73)*5/9); k
    36/(20*sqrt(73) + 36*I*sqrt(3) + 27)
    sage: N(k)
    0.165495678130644 - 0.0521492082074256*I
    sage: N(k,30)      # 30 "bits"
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{36}{20 \, \sqrt{73} + 36 i \, \sqrt{3} + 27}

Instalación
============

Si no tienes instalado Sage en tu computador y sólo quieres
probar algunos comandos, usa la versión en línea en http://sagecell.sagemath.org.

Mira la Guía De Instalación Para Sage en la sección de documentación de la
página web principal de [Sage]_ para obtener instrucciones sobre cómo instalar
Sage en tu computador. Aquí hacemos simplemente dos comentarios:


#. El archivo de descarga de Sage viene con "baterías incluidas". En otras
   palabras, aunque Sage utiliza Python, IPython, PARI, GAP, Singular,
   Maxima, NTL, GMP, etc., no necesitas instalarlos por separado
   pues ya están incluidos con la distribución de Sage.
   Sin embargo, para utilizar ciertas características de Sage, por ejemplo,
   Macaulay o KASH, debes tener los programas
   pertinentes ya instalados en tu computador.

#. La versión binaria precompilada de Sage (que se encuentra en el
   sitio web de Sage) puede ser más rápida y fácil de instalar que la
   versión en código fuente. Sólo desempaqueta el archivo y ejecuta ``sage``.


#. Si quieres utilizar el paquete SageTeX (el cual te permite insertar
   los resultados de tus cálculos con Sage en un archivo LaTeX),
   necesitarás hacerle conocer SageTeX a tu distribución de TeX.
   Para hacer esto, consulta la sección
   "Haciendo que TeX conozca a SageTeX" en la guía de instalación de Sage
   `Sage installation guide <http://doc.sagemath.org/html/en/installation/index.html>`_
   (`Este enlace
   <../../en/installation/index.html>`_ debería llevarte a tu copia
   local de la guía de instalación). Es bastante sencillo: sólo
   necesitas establecer una variable de entorno o copiar un solo archivo
   en un directorio en el que TeX va a buscar.

   La documentación para usar SageTeX se encuentra en
   ``$SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/``, donde
   "``$SAGE_ROOT``" se refiere al directorio donde Sage está instalado --
   por ejemplo, ``/opt/sage-9.6``.


Formas de usar Sage
===================

Puedes usar Sage de varias maneras.


-  **Interfaz gráfico del Notebook:** iniciar ``sage -n jupyter``; leer
   `Jupyter documentation on-line <https://jupyter-notebook.readthedocs.io/en/latest/notebook.html>`_,

-  **Línea de comandos interactiva:**,

-  **Programas:** Escribiendo programas compilados e interpretados en
   Sage y

-  **Scripts:** Escribiendo scripts (archivos de órdenes) independientes en Python
   que utilizan la biblioteca Sage.


Metas a largo plazo de Sage
===========================

-  **Útil**: La audiencia a la que está destinado Sage son los estudiantes de matemáticas
   (desde la secundaria hasta la universidad), profesores y matemáticos (para la investigación).
   El objetivo es proveer un software que pueda usarse para explorar y experimentar con construcciones
   matemáticas en álgebra, geometría, teoría de números, cálculo, computación numérica, etc.
   Sage facilita la experimentación interactiva con objetos matemáticos.

-  **Eficiente:** Queremos que sea rápido. Sage utiliza software maduro y altamente
   optimizado: GMP, PARI, GAP y NTL, por lo que es muy rápido en ciertas operaciones.

-  **Libre y de código abierto:** El código fuente debe ser legible y
   libremente disponible, de modo que los usuarios puedan entender qué está
   haciendo realmente el sistema y así poder extenderlo fácilmente. Tal como los matemáticos logran
   un entendimiento más profundo de un teorema al leerlo cuidadosamente o, por lo
   menos, al echarle una ojeada a la prueba, la gente que efectúa cálculos debe ser capaz de comprender
   cómo funcionan los cálculos leyendo el código fuente documentado.
   Si utilizas Sage para hacer cálculos en un artículo que vas a publicar,
   puedes estar seguro que tus lectores siempre tendrán libre acceso
   a Sage y a todo su código fuente, y hasta se te permite archivar y
   re-distribuir la versión de Sage que usaste.

-  **Fácil de compilar:** Sage tiene que ser fácil de compilar desde el
   código fuente para los usuarios de Linux, OS X y Windows. Esto provee
   a los usuarios de una mayor flexibilidad para que modifiquen el sistema.

-  **Cooperación con otros programas:** Sage debe proveer interfaces robustos a la mayoría de
   sistemas algebraicos de cómputo, incluyendo PARI, GAP, Singular, Maxima,
   KASH, Magma, Maple y Mathematica. Sage pretende unificar y extender
   el software matemático existente.

-  **Bien documentado:** Debemos proveer un tutorial, una guía de programación,
   un manual de referencia y documentos sobre cómo hacer cosas específicas,
   con numerosos ejemplos y discusiones de las bases matemáticas.

-  **Extensible:** Debe ser posible definir nuevos tipos de datos o derivar de
   tipos incorporados y utilizar código escrito en una amplia gama de lenguajes.

-  **Fácil de usar**: Debe de ser fácil comprender cual
   funcionalidad se ha provisto para un objeto dado y examinar
   la documentación y el código fuente, así como alcanzar un alto nivel
   de soporte al usuario.

.. [PyT] El Tutorial De Python,
   https://docs.python.org/es/3/tutorial/

.. [PyB] La Guía Para Principiantes De Python,
   https://wiki.python.org/moin/BeginnersGuide

.. [Sage] Sage, https://www.sagemath.org

