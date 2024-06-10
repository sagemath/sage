.. sage-doctest: needs sage.plot sage.symbolic

.. _section-plot:

Plotten
=======

Sage kann zwei- und dreidimensionale Plots erzeugen.

Zweidimensionale Plots
----------------------

Sage kann in zwei Dimensionen Kreise, Linien und Polygone zeichnen,
sowie Plots von Funktionen in kartesischen Koordinaten und Plots in
Polarkoordinaten, Konturplots und Plots von Vektorfeldern. Wir geben
davon im Folgenden einige Beispiele an. Für weitere Beispiele zum
Plotten mit Sage lesen Sie :ref:`section-systems` und
:ref:`section-maxima`, sowie die `Sage Constructions
<http://doc.sagemath.org/html/en/constructions/>`_ Dokumentation.

Dieser Befehl erstellt einen gelben Kreis vom Radius 1 mit dem
Ursprung als Zentrum:

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0))
    Graphics object consisting of 1 graphics primitive

Sie können auch einen ausgefüllten Kreis erzeugen:

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0), fill=True)
    Graphics object consisting of 1 graphics primitive

Sie können einen Kreis auch erstellen, indem Sie ihn einer Variable
zuweisen; so wird kein Plot gezeigt.

::

    sage: c = circle((0,0), 1, rgbcolor=(1,1,0))

Um den Plot zu zeigen, benutzen Sie ``c.show()`` oder ``show(c)`` wie
folgt:

.. link

::

    sage: c.show()

Alternativ führt das Auswerten von ``c.save('filename.png')`` dazu,
dass der Plot in der angegebenen Datei gespeichert wird.

Noch sehen diese 'Kreise' jedoch eher wie Ellipsen aus, da die Achsen
unterschiedlich skaliert sind. Sie können dies korrigieren:

.. link

::

    sage: c.show(aspect_ratio=1)

Der Befehl ``show(c, aspect_ratio=1)`` erreicht das Gleiche. Sie
können das Bild auch speichern, indem Sie ``c.save('filename.png',
aspect_ratio=1)`` benutzen.

Es ist einfach elementare Funktionen zu plotten:

::

    sage: plot(cos, (-5,5))
    Graphics object consisting of 1 graphics primitive

Sobald Sie einen Variablennamen angegeben haben, können Sie
parametrische Plots erzeugen:

::

    sage: x = var('x')
    sage: parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    Graphics object consisting of 1 graphics primitive

Es ist wichtig zu beachten, dass sich die Achsen eines Plots nur
schneiden, wenn sich der Ursprung im angezeigten Bildbereich des
Graphen befindet und ab einer bestimmten Größe der Werte wird die
wissenschaftliche Notation benutzt:
::

    sage: plot(x^2,(x,300,500))
    Graphics object consisting of 1 graphics primitive

Sie können mehrere Plots zusammenfügen indem Sie diese addieren:

::

    sage: x = var('x')
    sage: p1 = parametric_plot((cos(x),sin(x)),(x,0,2*pi),rgbcolor=hue(0.2))
    sage: p2 = parametric_plot((cos(x),sin(x)^2),(x,0,2*pi),rgbcolor=hue(0.4))
    sage: p3 = parametric_plot((cos(x),sin(x)^3),(x,0,2*pi),rgbcolor=hue(0.6))
    sage: show(p1+p2+p3, axes=false)

Eine gute Möglichkeit ausgefüllte Figuren zu erstellen ist, eine Liste
von Punkten zu erzeugen (``L`` im folgenden Beispiel) und dann den
``polygon`` Befehl zu verwenden um die Figur mit dem, durch die Punkte
bestimmten, Rand zu zeichnen. Zum Beispiel ist hier ein grünes Deltoid:

::

    sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),
    ....:     2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,3/4,1/2))
    sage: p
    Graphics object consisting of 1 graphics primitive

Geben Sie ``show(p, axes=false)`` ein, um dies ohne Achsen zu sehen.

Sie können auch Text zu einem Plot hinzufügen:

::

    sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),
    ....:     6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))
    sage: t = text("hypotrochoid", (5,4), rgbcolor=(1,0,0))
    sage: show(p+t)

Analysis Lehrer zeichnen häufig den folgenden Plot an die Tafel:
nicht nur einen Zweig von der arcsin Funktion, sondern mehrere, also den
Plot von  :math:`y=\sin(x)` für :math:`x` zwischen :math:`-2\pi` und
:math:`2\pi`, an der 45 Grad Linie gespiegelt. Der folgende Sage
Befehl erzeugt dies:

::

    sage: v = [(sin(x),x) for x in srange(-2*float(pi),2*float(pi),0.1)]
    sage: line(v)
    Graphics object consisting of 1 graphics primitive

Da die Tangensfunktion einen größeren Wertebereich als die
Sinusfunktion besitzt, sollten Sie, falls Sie den gleichen Trick
anwenden um die Inverse der Tangensfunktion zu plotten, das Minimum
und Maximum der Koordinaten für die *x*-Achse ändern:

::

    sage: v = [(tan(x),x) for x in srange(-2*float(pi),2*float(pi),0.01)]
    sage: show(line(v), xmin=-20, xmax=20)

Sage berechnet auch Plots in Polarkoordinaten, Konturplots, Plots von
Vektorfeldern (für besondere Arten von Funktionen). Hier ist ein
Beispiel eines Konturplots:

::

    sage: f = lambda x,y: cos(x*y)
    sage: contour_plot(f, (-4, 4), (-4, 4))
    Graphics object consisting of 1 graphics primitive

Dreidimensionale Plots
----------------------

Sage kann auch dazu verwendet werden dreidimensionale Plots zu zeichnen.
Sowohl im Notebook, als auch von der Kommandozeile aus werden diese
Plots standardmäßig mit den Open-Source-Paket [ThreeJS]_ angezeigt,
welches interaktives Drehen und Zoomen der Grafik mit Hilfe der
Maus unterstützt.

Benutzen Sie ``plot3d`` um eine Funktion der Form `f(x, y) = z` zu zeichnen:

::

    sage: x, y = var('x,y')
    sage: plot3d(x^2 + y^2, (x,-2,2), (y,-2,2))
    Graphics3d Object

Alternativ können Sie auch ``parametric_plot3d`` verwenden um eine
parametrisierte Fläche zu zeichnen, wobei jede der Variablen `x, y, z`
durch eine Funktion einer oder zweier Variablen bestimmt ist. (Die
Argumente sind typischerweise `u` und `v`). Der vorherige Plot kann
wie folgt parametrisiert angegeben werden:

::

    sage: u, v = var('u, v')
    sage: f_x(u, v) = u
    sage: f_y(u, v) = v
    sage: f_z(u, v) = u^2 + v^2
    sage: parametric_plot3d([f_x, f_y, f_z], (u, -2, 2), (v, -2, 2))
    Graphics3d Object

Die dritte Möglichkeit eine 3D Oberfläche zuplotten ist
``implicit_plot3d``, dies zeichnet eine Kontur einer Funktion mit
`f(x, y, z) = 0` (so wird eine Punktmenge definiert). Wir können die
Sphäre mithilfe einer klassischen Formel zeichnen:

::

    sage: x, y, z = var('x, y, z')
    sage: implicit_plot3d(x^2 + y^2 + z^2 - 4, (x,-2, 2), (y,-2, 2), (z,-2, 2))
    Graphics3d Object

Hier sind noch ein paar Beispiele:

`Whitneys Regenschirm <http://en.wikipedia.org/wiki/Whitney_umbrella>`__:

::

    sage: u, v = var('u,v')
    sage: fx = u*v
    sage: fy = u
    sage: fz = v^2
    sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -1, 1),
    ....:   frame=False, color="yellow")
    Graphics3d Object

Die `Kreuz-Kappe <http://de.wikipedia.org/wiki/Kreuzhaube>`__:

::

    sage: u, v = var('u,v')
    sage: fx = (1+cos(v))*cos(u)
    sage: fy = (1+cos(v))*sin(u)
    sage: fz = -tanh((2/3)*(u-pi))*sin(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....:   frame=False, color="red")
    Graphics3d Object

Ein gedrehter Torus:

::

    sage: u, v = var('u,v')
    sage: fx = (3+sin(v)+cos(u))*cos(2*v)
    sage: fy = (3+sin(v)+cos(u))*sin(2*v)
    sage: fz = sin(u)+2*cos(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ....:   frame=False, color="red")
    Graphics3d Object

Die `Lemniskate <http://de.wikipedia.org/wiki/Lemniskate>`__:

::

    sage: x, y, z = var('x,y,z')
    sage: f(x, y, z) = 4*x^2 * (x^2 + y^2 + z^2 + z) + y^2 * (y^2 + z^2 - 1)
    sage: implicit_plot3d(f, (x, -0.5, 0.5), (y, -1, 1), (z, -1, 1))
    Graphics3d Object
