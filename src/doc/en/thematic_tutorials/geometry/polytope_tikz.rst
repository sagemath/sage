.. -*- coding: utf-8 -*-

.. linkall

.. _polytikz:

Draw polytopes in LaTeX using TikZ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>


It is sometimes very helpful to draw 3-dimensional polytopes in a
paper. TikZ is a very versatile tool to draw in scientific documents
and Sage can deal easily with 3-dimensional polytopes. Finally sagetex
makes everything work together nicely between Sage, TikZ and
LaTeX. Since version 6.3 of Sage, there is a function for (projection
of) polytopes to output a TikZ picture of the polytope. Since version 9.8 of
SageMath, the tikz output can be a ``TikzPicture`` object from the sage module
``sage.misc.latex_standalone``. This short tutorial shows how it all works.

Instructions
""""""""""""

To put an image of a 3D-polytope in LaTeX using TikZ and Sage, simply follow the instructions:

- Install `SageTex <http://doc.sagemath.org/html/en/tutorial/sagetex.html>`_ (optional but recommended!)
- Put ``\usepackage{tikz}`` in the preamble of your article
- Open Sage and change the directory to your article's by the command ``cd /path/to/article``
- Input your polytope, called P for example, to Sage
- Visualize the polytope P using the command ``P.show(aspect_ratio=1)``
- This will open an interactive view in your default browser, where you can rotate the polytope.
- Once the desired view angle is found, click on the information icon in the lower right-hand corner and select *Get Viewpoint*. This will copy a string of the form '[x,y,z],angle' to your local clipboard.
- Go back to Sage and type ``Img = P.tikz([x,y,z],angle,output_type='LatexExpr')``. You can paste the string here to save some typing.
- *Img* now contains a Sage object of type ``LatexExpr`` containing the raw TikZ picture of your polytope.

Alternatively, you can save the tikz image to a file, by doing

.. CODE-BLOCK:: python

  Img = P.tikz([x,y,z], angle, output_type='TikzPicture')
  Img.tex('Img_poly.tex')
  Img.tex('Img_poly.tex', content_only=True)
  Img.pdf('Img_poly.pdf')

.. end of output

Then in the pwd (present working directory of sage, the one of your article)
you will have two files named ``Img_poly.tex`` and ``Img_poly.pdf`` containing the
tikzpicture of the polytope ``P``.

Customization
"""""""""""""

You can customize the polytope using the following options in the command ``P.tikz()``

- ``scale`` : positive number to scale the polytope,
- ``edge_color`` : string (default: ``blue!95!black``) representing colors which tikz recognize,
- ``facet_color`` : string (default: ``blue!95!black``) representing colors which tikz recognize,
- ``vertex_color`` : string (default: ``green``) representing colors which tikz recognize,
- ``opacity`` : real number (default: ``0.8``) between 0 and 1 giving the opacity of the front facets,
- ``axis`` : Boolean (default: ``False``) draw the axes at the origin or not.
- ``output_type`` : string (default: ``None``) ``None``, ``'LatexExpr'`` or
  ``'TikzPicture'``, the type of the output. Since SageMath 9.8, the value ``None`` is deprecated
  as the default value will soon be changed from ``'LatexExpr'`` to ``'TikzPicture'``.

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

    Img = P.tikz([674,108,-731], 112, output_type='TikzPicture')

.. end of output

Note: the ``output_type='TikzPicture'`` is necessary since SagMath 9.8 to avoid
a deprecation warning message since the default output type will soon change
from a ``LatexExpr`` (Python str) to a ``TikzPicture`` object (allowing more
versatility, like being able to view it directly in the Jupyter notebook).

Or you may want to customize using the command

::

    Img = P.tikz([674,108,-731],112,scale=2, edge_color='orange',facet_color='red',vertex_color='blue',opacity=0.4, output_type='TikzPicture')

.. end of output

Further, you may want to edit deeper the style of the polytope, directly inside the tikzpicture. For example, line 6-9 in the tikzpicture reads:

.. CODE-BLOCK:: latex

  back/.style={loosely dotted, thin},
  edge/.style={color=orange, thick},
  facet/.style={fill=red,fill opacity=0.400000},
  vertex/.style={inner sep=1pt,circle,draw=blue!25!black,fill=blue!75!black,thick,anchor=base}]

.. end of output


It is also possible to replace it by the following 4 lines (and adding ``\usetikzlibrary{shapes}`` in the preamble)

.. CODE-BLOCK:: latex

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

.. CODE-BLOCK:: latex

  \sagestr{(polytopes.permutahedron(4)).tikz([4,5,6],45,scale=0.75, facet_color='red',vertex_color='yellow',opacity=0.3, output_type='LatexExpr')}

.. end of output

2) You may create the following tex commands

.. CODE-BLOCK:: latex

  \newcommand{\polytopeimg}[4]{\sagestr{(#1).tikz(#2,#3,#4,output_type='LatexExpr')}}
  \newcommand{\polytopeimgopt}[9]{\sagestr{(#1).tikz(#2,#3,#4,#5,#6,#7,#8,#9,output_type='LatexExpr')}}

.. end of output

in your preamble and use them with a sagesilent in your article:

.. CODE-BLOCK:: latex

  \begin{sagesilent}
  Polytope = polytopes.great_rhombicuboctahedron()
  \end{sagesilent}

.. end of output

.. CODE-BLOCK:: latex

  \polytopeimg{Polytope}{[276,-607,-746]}{102}{1}
  \polytopeimgopt{Polytope}{view=[-907,379,183]}{angle=129}{scale=2}{edge_color='red'}{facet_color='yellow'}{vertex_color='blue'}{opacity=0.6}{axis=False}

.. end of output

Then, run pdflatex, execute Sage on the file ``article_name.sagetex.sage`` and run pdflatex again.

For more information on SageTeX, see the tutorial http://doc.sagemath.org/html/en/tutorial/sagetex.html.

