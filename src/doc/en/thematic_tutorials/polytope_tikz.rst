.. -*- coding: utf-8 -*-

.. linkall

.. _polytikz:

Draw polytopes in LateX using TikZ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.huji.ac.il>


It is sometimes very helpful to draw 3-dimensional polytopes in a
paper. TikZ is a very versatile tool to draw in scientific documents
and Sage can deal easily with 3-dimensional polytopes. Finally sagetex
makes everything work together nicely between Sage, TikZ and
LaTeX. Since version 6.3 of Sage, there is a function for (projection
of) polytopes to output a TikZ picture of the polytope.  This short
tutorial shows how it all works.

Instructions
""""""""""""

To put an image of a 3d-polytope in LaTeX using TikZ and Sage, simply follow the instructions:

- Install `SageTex <http://doc.sagemath.org/html/en/tutorial/sagetex.html>`_ (optional but recommended!)
- Put ``\usepackage{tikz}`` in the preamble of your article
- Open Sage and change the directory to your article's by the command ``cd /path/to/article``
- Input your polytope, called P for example, to Sage
- Visualize the polytope P using the command ``P.show(aspect_ratio=1)``
- This will open an interactive viewer named Jmol, in which you can rotate the polytope. Once the wished view angle is found, right click on the image and select *Console*
- In the dialog box click the button *State*
- Scroll up to the line starting with *moveto*
- It reads something like ``moveto 0.0 {x y z angle} scale``
- Go back to Sage and type ``Img = P.projection().tikz([x,y,z],angle)``
- *Img* now contains a Sage object of type ``LatexExpr`` containing the raw TikZ picture of your polytope

Then, you can either copy-paste it to your article by typing ``Img`` in Sage or save it to a file, by doing

::

  f = open('Img_poly.tex','w')
  f.write(Img)
  f.close()

.. end of output

Then in the pwd (present working directory of sage, the one of your article) 
you will have a file named ``Img_poly.tex`` containing the tikzpicture of your polytope.

Customization
"""""""""""""

You can customize the polytope using the following options in the command ``P.tikz()``

- ``scale`` : positive number to scale the polytope,
- ``edge_color`` : string (default: ``blue!95!black``) representing colors which tikz recognize,
- ``facet_color`` : string (default: ``blue!95!black``) representing colors which tikz recognize,
- ``vertex_color`` : string (default: ``green``) representing colors which tikz recognize,
- ``opacity`` : real number (default: ``0.8``) between 0 and 1 giving the opacity of the front facets,
- ``axis`` : Boolean (default: ``False``) draw the axes at the origin or not.

Examples
""""""""

Let's say you want to draw the polar dual of the following (nice!) polytope given by the following list of vertices:

``[[1,0,1],[1,0,0],[1,1,0],[0,0,-1],[0,1,0],[-1,0,0],[0,1,1],[0,0,1],[0,-1,0]]``

In Sage, you type:

::

    P = Polyhedron(vertices=[[1,0,1],[1,0,0],[1,1,0],[0,0,-1],[0,1,0],[-1,0,0],[0,1,1],[0,0,1],[0,-1,0]]).polar()

.. end of output

Then, you visualize the polytope by typing ``P.show(aspect_ratio=1)``

When you found a good angle, follow the above procedure to obtain the values 
[674,108,-731] and angle=112, for example.

::

    Img = P.projection().tikz([674,108,-731],112)

.. end of output

Or you may want to customize using the command

::

    Img = P.projection().tikz([674,108,-731],112,scale=2, edge_color='orange',facet_color='red',vertex_color='blue',opacity=0.4)

.. end of output

Further, you may want to edit deeper the style of the polytope, directly inside the tikzpicture. For example, line 6-9 in the tikzpicture reads:

::

  back/.style={loosely dotted, thin},
  edge/.style={color=orange, thick},
  facet/.style={fill=red,fill opacity=0.400000},
  vertex/.style={inner sep=1pt,circle,draw=blue!25!black,fill=blue!75!black,thick,anchor=base}]

.. end of output


It is also possible to replace it by the following 4 lines (and adding ``\usetikzlibrary{shapes}`` in the preamble)

::

  back/.style={loosely dashed,line width=2pt},
  edge/.style={color=yellow, line width=2pt},
  facet/.style={fill=cyan,fill opacity=0.400000},
  vertex/.style={inner sep=4pt,star,star points=7,draw=blue!75!white,fill=blue!85!white,thick,anchor=base}]

.. end of output

Finally, you may want to tweak your picture my adding labels, elements on 
vertices, edges, facets, etc.

Automatize using SageTex
""""""""""""""""""""""""

For this you need to put

``\usepackage{sagetex}``

in the preamble of your article

There are different ways to use sagetex and you may create your own. Here are 
some possibilities.

1) You can directly type in a sagestr in the article:

::

  \sagestr{(polytopes.permutahedron(4)).projection().tikz([4,5,6],45,scale=0.75, facet_color='red',vertex_color='yellow',opacity=0.3)}

.. end of output

2) You may create the following tex commands

::

  \newcommand{\polytopeimg}[4]{\sagestr{(#1).projection().tikz(#2,#3,#4)}}
  \newcommand{\polytopeimgopt}[9]{\sagestr{(#1).projection().tikz(#2,#3,#4,#5,#6,#7,#8,#9)}}

.. end of output

in your preamble and use them with a sagesilent in your article:

::

  \begin{sagesilent}
  Polytope = polytopes.great_rhombicuboctahedron()
  \end{sagesilent}

.. end of output

::

  \polytopeimg{Polytope}{[276,-607,-746]}{102}{1}
  \polytopeimgopt{Polytope}{view=[-907,379,183]}{angle=129}{scale=2}{edge_color='red'}{facet_color='yellow'}{vertex_color='blue'}{opacity=0.6}{axis=False}

.. end of output

Then, run pdflatex, execute Sage on the file ``article_name.sagetex.sage`` and run pdflatex again.

For more information on SageTeX, see the tutorial http://doc.sagemath.org/html/en/tutorial/sagetex.html.

