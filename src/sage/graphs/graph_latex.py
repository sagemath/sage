# sage.doctest: needs sage.plot
r"""
LaTeX options for graphs

This module provides a class to hold, manipulate and employ various options for
rendering a graph in LaTeX, in addition to providing the code that actually
generates a LaTeX representation of a (combinatorial) graph.

AUTHORS:

- Rob Beezer (2009-05-20): :class:`~sage.graphs.graph_latex.GraphLatex` class
- Fidel Barerra Cruz (2009-05-20): ``tkz-graph`` commands to render a graph
- Nicolas M. Thiéry (2010-02): dot2tex/graphviz interface
- Rob Beezer (2010-05-29): Extended range of ``tkz-graph`` options

LaTeX Versions of Graphs
-------------------------------------

.. image:: ../../media/heawood-graph-latex.png
   :align: center

Many mathematical objects in Sage have LaTeX representations, and graphs are no
exception.  For a graph ``g``, the command ``view(g)``, issued at the Sage
command line or in the notebook, will create a graphic version of ``g``.
Similarly, ``latex(g)`` will return a (long) string that is a representation of
the graph in LaTeX.  Other ways of employing LaTeX in Sage, such as ``%latex``
in a notebook cell, or the Typeset checkbox in the notebook, will handle ``g``
appropriately.

Support through the ``tkz-graph`` package is by Alain Matthes, the author of
``tkz-graph``, whose work can be found at https://www.ctan.org/pkg/tkz-graph.

The range of possible options for customizing the appearance of a graph are
carefully documented at :meth:`sage.graphs.graph_latex.GraphLatex.set_option`.
As a broad overview, the following options are supported:

- Pre-built Styles: the pre-built styles of the tkz-graph package provide nice
  drawings quickly
- Dimensions: can be specified in natural units, then uniformly scaled after
  design work
- Vertex Colors: the perimeter and fill color for vertices can be specified,
  including on a per-vertex basis
- Vertex Shapes: may be circles, shaded spheres, rectangles or diamonds,
  including on a per-vertex basis
- Vertex Sizes: may be specified as minimums, and will automatically sized to
  contain vertex labels, including on a per-vertex basis
- Vertex Labels: can use latex formatting, and may have their colors specified,
  including on a per-vertex basis
- Vertex Label Placement: can be interior to the vertex, or external at a
  configurable location
- Edge Colors: a solid color with or without a second color down the middle, on
  a per-edge basis
- Edge Thickness: can be set, including on a per-edge basis
- Edge Labels: can use latex formatting, and may have their colors specified,
  including on a per-edge basis
- Edge Label Placement: can be to the left, right, above, below, inline, and
  then sloped or horizontal
- Digraph Edges: are slightly curved, with arrowheads
- Loops: may be specified by their size, and with a direction equaling one of
  the four compass points

To use LaTeX in Sage you of course need a working TeX installation and it will
work best if you have the ``dvipng`` and ``magick`` utilities.  For graphs you
need the ``tkz-graph.sty`` and ``tkz-berge.sty`` style files of the tkz-graph
package.  TeX, dvipng, and convert should be widely available through package
managers or installers.  You may need to install the tkz-graph style files in
the appropriate locations, a task beyond the scope of this introduction.
Primary locations for these programs are:

- TeX: http://ctan.org/
- dvipng: http://sourceforge.net/projects/dvipng/
- magick: http://www.imagemagick.org (the ImageMagick suite)
- tkz-graph: https://www.ctan.org/pkg/tkz-graph

Customizing the output is accomplished in several ways.  Suppose ``g`` is a
graph, then ``g.set_latex_options()`` can be used to efficiently set or modify
various options.  Setting individual options, or querying options, can be
accomplished by first using a command like ``opts = g.latex_options()`` to
obtain a :class:`sage.graphs.graph_latex.GraphLatex` object which has several
methods to set and retrieve options.

Here is a minimal session demonstrating how to use these features. The following
setup should work in the notebook or at the command-line.::

    sage: H = graphs.HeawoodGraph()
    sage: H.set_latex_options(
    ....:     graphic_size=(5,5),
    ....:     vertex_size=0.2,
    ....:     edge_thickness=0.04,
    ....:     edge_color='green',
    ....:     vertex_color='green',
    ....:     vertex_label_color='red'
    ....: )

At this point, ``view(H)`` should call ``pdflatex`` to process the string
created by ``latex(H)`` and then display the resulting graphic.

To use this image in a LaTeX document, you could of course just copy and save
the resulting graphic.  However, the ``latex()`` command will produce the
underlying LaTeX code, which can be incorporated into a standalone LaTeX
document.::

    sage: from sage.graphs.graph_latex import check_tkz_graph
    sage: check_tkz_graph()  # random - depends on TeX installation
    sage: latex(H)
    \begin{tikzpicture}
    \definecolor{cv0}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv0}{rgb}{1.0,1.0,1.0}
    \definecolor{clv0}{rgb}{1.0,0.0,0.0}
    \definecolor{cv1}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv1}{rgb}{1.0,1.0,1.0}
    \definecolor{clv1}{rgb}{1.0,0.0,0.0}
    \definecolor{cv2}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv2}{rgb}{1.0,1.0,1.0}
    \definecolor{clv2}{rgb}{1.0,0.0,0.0}
    \definecolor{cv3}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv3}{rgb}{1.0,1.0,1.0}
    \definecolor{clv3}{rgb}{1.0,0.0,0.0}
    \definecolor{cv4}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv4}{rgb}{1.0,1.0,1.0}
    \definecolor{clv4}{rgb}{1.0,0.0,0.0}
    \definecolor{cv5}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv5}{rgb}{1.0,1.0,1.0}
    \definecolor{clv5}{rgb}{1.0,0.0,0.0}
    \definecolor{cv6}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv6}{rgb}{1.0,1.0,1.0}
    \definecolor{clv6}{rgb}{1.0,0.0,0.0}
    \definecolor{cv7}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv7}{rgb}{1.0,1.0,1.0}
    \definecolor{clv7}{rgb}{1.0,0.0,0.0}
    \definecolor{cv8}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv8}{rgb}{1.0,1.0,1.0}
    \definecolor{clv8}{rgb}{1.0,0.0,0.0}
    \definecolor{cv9}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv9}{rgb}{1.0,1.0,1.0}
    \definecolor{clv9}{rgb}{1.0,0.0,0.0}
    \definecolor{cv10}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv10}{rgb}{1.0,1.0,1.0}
    \definecolor{clv10}{rgb}{1.0,0.0,0.0}
    \definecolor{cv11}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv11}{rgb}{1.0,1.0,1.0}
    \definecolor{clv11}{rgb}{1.0,0.0,0.0}
    \definecolor{cv12}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv12}{rgb}{1.0,1.0,1.0}
    \definecolor{clv12}{rgb}{1.0,0.0,0.0}
    \definecolor{cv13}{rgb}{0.0,0.502,0.0}
    \definecolor{cfv13}{rgb}{1.0,1.0,1.0}
    \definecolor{clv13}{rgb}{1.0,0.0,0.0}
    \definecolor{cv0v1}{rgb}{0.0,0.502,0.0}
    \definecolor{cv0v5}{rgb}{0.0,0.502,0.0}
    \definecolor{cv0v13}{rgb}{0.0,0.502,0.0}
    \definecolor{cv1v2}{rgb}{0.0,0.502,0.0}
    \definecolor{cv1v10}{rgb}{0.0,0.502,0.0}
    \definecolor{cv2v3}{rgb}{0.0,0.502,0.0}
    \definecolor{cv2v7}{rgb}{0.0,0.502,0.0}
    \definecolor{cv3v4}{rgb}{0.0,0.502,0.0}
    \definecolor{cv3v12}{rgb}{0.0,0.502,0.0}
    \definecolor{cv4v5}{rgb}{0.0,0.502,0.0}
    \definecolor{cv4v9}{rgb}{0.0,0.502,0.0}
    \definecolor{cv5v6}{rgb}{0.0,0.502,0.0}
    \definecolor{cv6v7}{rgb}{0.0,0.502,0.0}
    \definecolor{cv6v11}{rgb}{0.0,0.502,0.0}
    \definecolor{cv7v8}{rgb}{0.0,0.502,0.0}
    \definecolor{cv8v9}{rgb}{0.0,0.502,0.0}
    \definecolor{cv8v13}{rgb}{0.0,0.502,0.0}
    \definecolor{cv9v10}{rgb}{0.0,0.502,0.0}
    \definecolor{cv10v11}{rgb}{0.0,0.502,0.0}
    \definecolor{cv11v12}{rgb}{0.0,0.502,0.0}
    \definecolor{cv12v13}{rgb}{0.0,0.502,0.0}
    %
    \Vertex[style={minimum size=0.2cm,draw=cv0,fill=cfv0,text=clv0,shape=circle},LabelOut=false,L=\hbox{$0$},x=2.5cm,y=5.0cm]{v0}
    \Vertex[style={minimum size=0.2cm,draw=cv1,fill=cfv1,text=clv1,shape=circle},LabelOut=false,L=\hbox{$1$},x=1.3874cm,y=4.7524cm]{v1}
    \Vertex[style={minimum size=0.2cm,draw=cv2,fill=cfv2,text=clv2,shape=circle},LabelOut=false,L=\hbox{$2$},x=0.4952cm,y=4.0587cm]{v2}
    \Vertex[style={minimum size=0.2cm,draw=cv3,fill=cfv3,text=clv3,shape=circle},LabelOut=false,L=\hbox{$3$},x=0.0cm,y=3.0563cm]{v3}
    \Vertex[style={minimum size=0.2cm,draw=cv4,fill=cfv4,text=clv4,shape=circle},LabelOut=false,L=\hbox{$4$},x=0.0cm,y=1.9437cm]{v4}
    \Vertex[style={minimum size=0.2cm,draw=cv5,fill=cfv5,text=clv5,shape=circle},LabelOut=false,L=\hbox{$5$},x=0.4952cm,y=0.9413cm]{v5}
    \Vertex[style={minimum size=0.2cm,draw=cv6,fill=cfv6,text=clv6,shape=circle},LabelOut=false,L=\hbox{$6$},x=1.3874cm,y=0.2476cm]{v6}
    \Vertex[style={minimum size=0.2cm,draw=cv7,fill=cfv7,text=clv7,shape=circle},LabelOut=false,L=\hbox{$7$},x=2.5cm,y=0.0cm]{v7}
    \Vertex[style={minimum size=0.2cm,draw=cv8,fill=cfv8,text=clv8,shape=circle},LabelOut=false,L=\hbox{$8$},x=3.6126cm,y=0.2476cm]{v8}
    \Vertex[style={minimum size=0.2cm,draw=cv9,fill=cfv9,text=clv9,shape=circle},LabelOut=false,L=\hbox{$9$},x=4.5048cm,y=0.9413cm]{v9}
    \Vertex[style={minimum size=0.2cm,draw=cv10,fill=cfv10,text=clv10,shape=circle},LabelOut=false,L=\hbox{$10$},x=5.0cm,y=1.9437cm]{v10}
    \Vertex[style={minimum size=0.2cm,draw=cv11,fill=cfv11,text=clv11,shape=circle},LabelOut=false,L=\hbox{$11$},x=5.0cm,y=3.0563cm]{v11}
    \Vertex[style={minimum size=0.2cm,draw=cv12,fill=cfv12,text=clv12,shape=circle},LabelOut=false,L=\hbox{$12$},x=4.5048cm,y=4.0587cm]{v12}
    \Vertex[style={minimum size=0.2cm,draw=cv13,fill=cfv13,text=clv13,shape=circle},LabelOut=false,L=\hbox{$13$},x=3.6126cm,y=4.7524cm]{v13}
    %
    \Edge[lw=0.04cm,style={color=cv0v1,},](v0)(v1)
    \Edge[lw=0.04cm,style={color=cv0v5,},](v0)(v5)
    \Edge[lw=0.04cm,style={color=cv0v13,},](v0)(v13)
    \Edge[lw=0.04cm,style={color=cv1v2,},](v1)(v2)
    \Edge[lw=0.04cm,style={color=cv1v10,},](v1)(v10)
    \Edge[lw=0.04cm,style={color=cv2v3,},](v2)(v3)
    \Edge[lw=0.04cm,style={color=cv2v7,},](v2)(v7)
    \Edge[lw=0.04cm,style={color=cv3v4,},](v3)(v4)
    \Edge[lw=0.04cm,style={color=cv3v12,},](v3)(v12)
    \Edge[lw=0.04cm,style={color=cv4v5,},](v4)(v5)
    \Edge[lw=0.04cm,style={color=cv4v9,},](v4)(v9)
    \Edge[lw=0.04cm,style={color=cv5v6,},](v5)(v6)
    \Edge[lw=0.04cm,style={color=cv6v7,},](v6)(v7)
    \Edge[lw=0.04cm,style={color=cv6v11,},](v6)(v11)
    \Edge[lw=0.04cm,style={color=cv7v8,},](v7)(v8)
    \Edge[lw=0.04cm,style={color=cv8v9,},](v8)(v9)
    \Edge[lw=0.04cm,style={color=cv8v13,},](v8)(v13)
    \Edge[lw=0.04cm,style={color=cv9v10,},](v9)(v10)
    \Edge[lw=0.04cm,style={color=cv10v11,},](v10)(v11)
    \Edge[lw=0.04cm,style={color=cv11v12,},](v11)(v12)
    \Edge[lw=0.04cm,style={color=cv12v13,},](v12)(v13)
    %
    \end{tikzpicture}

EXAMPLES:

This example illustrates switching between the built-in styles when using the
tkz_graph format.::

    sage: g = graphs.PetersenGraph()
    sage: g.set_latex_options(tkz_style='Classic')
    sage: from sage.graphs.graph_latex import check_tkz_graph
    sage: check_tkz_graph()  # random - depends on TeX installation
    sage: latex(g)
    \begin{tikzpicture}
    \GraphInit[vstyle=Classic]
    ...
    \end{tikzpicture}
    sage: opts = g.latex_options()
    sage: opts
    LaTeX options for Petersen graph: {'tkz_style': 'Classic'}
    sage: g.set_latex_options(tkz_style = 'Art')
    sage: opts.get_option('tkz_style')
    'Art'
    sage: opts
    LaTeX options for Petersen graph: {'tkz_style': 'Art'}
    sage: latex(g)
    \begin{tikzpicture}
    \GraphInit[vstyle=Art]
    ...
    \end{tikzpicture}

This example illustrates using the optional dot2tex module::

    sage: g = graphs.PetersenGraph()
    sage: g.set_latex_options(format='dot2tex', prog='neato')
    sage: from sage.graphs.graph_latex import check_tkz_graph
    sage: check_tkz_graph()  # random - depends on TeX installation
    sage: latex(g)  # optional - dot2tex graphviz
    \begin{tikzpicture}[>=latex,line join=bevel,]
    ...
    \end{tikzpicture}

Among other things, this supports the flexible ``edge_options`` option
(see :meth:`sage.graphs.generic_graph.GenericGraph.graphviz_string`);
here we color in red all edges touching the vertex ``0``::

    sage: g = graphs.PetersenGraph()
    sage: g.set_latex_options(format='dot2tex', edge_options=lambda u_v_label: {"color": "red"} if u_v_label[0] == 0 else {})
    sage: latex(g)  # optional - dot2tex graphviz
    \begin{tikzpicture}[>=latex,line join=bevel,]
    ...
    \end{tikzpicture}


TESTS:

This graph will look horrible, but it illustrates (and tests) a great variety of
the possible options available through Sage's interface to the ``tkz-graph``
package.  So it is worth viewing this in the notebook to see the effects of
various defaults and choices.::

    sage: # needs sage.symbolic
    sage: var('x y u w')
    (x, y, u, w)
    sage: G = Graph(loops=True)
    sage: for i in range(5):
    ....:    for j in range(i+1, 5):
    ....:         G.add_edge((i, j), label=(x^i*y^j).expand())
    sage: G.add_edge((0,0), label=sin(u))
    sage: G.add_edge((4,4), label=w^5)
    sage: G.set_pos(G.layout_circular())
    sage: G.set_latex_options(
    ....:     units='in',
    ....:     graphic_size=(8,8),
    ....:     margins=(1,2,2,1),
    ....:     scale=0.5,
    ....:     vertex_color='0.8',
    ....:     vertex_colors={1:'aqua', 3:'y', 4:'#0000FF'},
    ....:     vertex_fill_color='blue',
    ....:     vertex_fill_colors={1:'green', 3:'b', 4:'#FF00FF'},
    ....:     vertex_label_color='brown',
    ....:     vertex_label_colors={0:'g',1:'purple',2:'#007F00'},
    ....:     vertex_shape='diamond',
    ....:     vertex_shapes={1:'rectangle', 2:'sphere', 3:'sphere', 4:'circle'},
    ....:     vertex_size=0.3,
    ....:     vertex_sizes={0:1.0, 2:0.3, 4:1.0},
    ....:     vertex_label_placements = {2:(0.6, 180), 4:(0,45)},
    ....:     edge_color='purple',
    ....:     edge_colors={(0,2):'g',(3,4):'red'},
    ....:     edge_fills=True,
    ....:     edge_fill_color='green',
    ....:     edge_label_colors={(2,3):'y',(0,4):'blue'},
    ....:     edge_thickness=0.05,
    ....:     edge_thicknesses={(3,4):0.2, (0,4):0.02},
    ....:     edge_labels=True,
    ....:     edge_label_sloped=True,
    ....:     edge_label_slopes={(0,3):False, (2,4):False},
    ....:     edge_label_placement=0.50,
    ....:     edge_label_placements={(0,4):'above', (2,3):'left', (0,0):'above', (4,4):'below'},
    ....:     loop_placement=(2.0, 'NO'),
    ....:     loop_placements={4:(8.0, 'EA')}
    ....: )
    sage: from sage.graphs.graph_latex import check_tkz_graph
    sage: check_tkz_graph()  # random - depends on TeX installation
    sage: print(latex(G))
    \begin{tikzpicture}
    \definecolor{cv0}{rgb}{0.8,0.8,0.8}
    \definecolor{cfv0}{rgb}{0.0,0.0,1.0}
    \definecolor{clv0}{rgb}{0.0,0.5,0.0}
    \definecolor{cv1}{rgb}{0.0,1.0,1.0}
    \definecolor{cfv1}{rgb}{0.0,0.502,0.0}
    \definecolor{clv1}{rgb}{0.502,0.0,0.502}
    \definecolor{cv2}{rgb}{0.8,0.8,0.8}
    \definecolor{cfv2}{rgb}{0.0,0.0,1.0}
    \definecolor{clv2}{rgb}{0.0,0.498,0.0}
    \definecolor{cv3}{rgb}{0.75,0.75,0.0}
    \definecolor{cfv3}{rgb}{0.0,0.0,1.0}
    \definecolor{clv3}{rgb}{0.6471,0.1647,0.1647}
    \definecolor{cv4}{rgb}{0.0,0.0,1.0}
    \definecolor{cfv4}{rgb}{1.0,0.0,1.0}
    \definecolor{clv4}{rgb}{0.6471,0.1647,0.1647}
    \definecolor{cv0v0}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv0v0}{rgb}{0.0,0.502,0.0}
    \definecolor{clv0v0}{rgb}{0.0,0.0,0.0}
    \definecolor{cv0v1}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv0v1}{rgb}{0.0,0.502,0.0}
    \definecolor{clv0v1}{rgb}{0.0,0.0,0.0}
    \definecolor{cv0v2}{rgb}{0.0,0.5,0.0}
    \definecolor{cfv0v2}{rgb}{0.0,0.502,0.0}
    \definecolor{clv0v2}{rgb}{0.0,0.0,0.0}
    \definecolor{cv0v3}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv0v3}{rgb}{0.0,0.502,0.0}
    \definecolor{clv0v3}{rgb}{0.0,0.0,0.0}
    \definecolor{cv0v4}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv0v4}{rgb}{0.0,0.502,0.0}
    \definecolor{clv0v4}{rgb}{0.0,0.0,1.0}
    \definecolor{cv1v2}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv1v2}{rgb}{0.0,0.502,0.0}
    \definecolor{clv1v2}{rgb}{0.0,0.0,0.0}
    \definecolor{cv1v3}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv1v3}{rgb}{0.0,0.502,0.0}
    \definecolor{clv1v3}{rgb}{0.0,0.0,0.0}
    \definecolor{cv1v4}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv1v4}{rgb}{0.0,0.502,0.0}
    \definecolor{clv1v4}{rgb}{0.0,0.0,0.0}
    \definecolor{cv2v3}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv2v3}{rgb}{0.0,0.502,0.0}
    \definecolor{clv2v3}{rgb}{0.75,0.75,0.0}
    \definecolor{cv2v4}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv2v4}{rgb}{0.0,0.502,0.0}
    \definecolor{clv2v4}{rgb}{0.0,0.0,0.0}
    \definecolor{cv3v4}{rgb}{1.0,0.0,0.0}
    \definecolor{cfv3v4}{rgb}{0.0,0.502,0.0}
    \definecolor{clv3v4}{rgb}{0.0,0.0,0.0}
    \definecolor{cv4v4}{rgb}{0.502,0.0,0.502}
    \definecolor{cfv4v4}{rgb}{0.0,0.502,0.0}
    \definecolor{clv4v4}{rgb}{0.0,0.0,0.0}
    %
    \Vertex[style={minimum size=0.5in,draw=cv0,fill=cfv0,text=clv0,shape=diamond},LabelOut=false,L=\hbox{$0$},x=1.75in,y=3.0in]{v0}
    \Vertex[style={minimum size=0.15in,draw=cv1,fill=cfv1,text=clv1,shape=rectangle},LabelOut=false,L=\hbox{$1$},x=0.5in,y=2.0451in]{v1}
    \Vertex[style={minimum size=0.15in,draw=cv2,fill=cfv2,text=clv2,shape=circle,shading=ball,line width=0pt,ball color=cv2,},LabelOut=true,Ldist=0.3in,Lpos=180.0,L=\hbox{$2$},x=0.9775in,y=0.5in]{v2}
    \Vertex[style={minimum size=0.15in,draw=cv3,fill=cfv3,text=clv3,shape=circle,shading=ball,line width=0pt,ball color=cv3,},LabelOut=false,L=\hbox{$3$},x=2.5225in,y=0.5in]{v3}
    \Vertex[style={minimum size=0.5in,draw=cv4,fill=cfv4,text=clv4,shape=circle},LabelOut=true,Ldist=0.0in,Lpos=45.0,L=\hbox{$4$},x=3.0in,y=2.0451in]{v4}
    %
    \Loop[dist=1.0in,dir=NO,style={color=cv0v0,double=cfv0v0},labelstyle={sloped,above,text=clv0v0,},label=\hbox{$\sin\left(u\right)$},](v0)
    \Edge[lw=0.025in,style={color=cv0v1,double=cfv0v1},labelstyle={sloped,pos=0.5,text=clv0v1,},label=\hbox{$y$},](v0)(v1)
    \Edge[lw=0.025in,style={color=cv0v2,double=cfv0v2},labelstyle={sloped,pos=0.5,text=clv0v2,},label=\hbox{$y^{2}$},](v0)(v2)
    \Edge[lw=0.025in,style={color=cv0v3,double=cfv0v3},labelstyle={pos=0.5,text=clv0v3,},label=\hbox{$y^{3}$},](v0)(v3)
    \Edge[lw=0.01in,style={color=cv0v4,double=cfv0v4},labelstyle={sloped,above,text=clv0v4,},label=\hbox{$y^{4}$},](v0)(v4)
    \Edge[lw=0.025in,style={color=cv1v2,double=cfv1v2},labelstyle={sloped,pos=0.5,text=clv1v2,},label=\hbox{$x y^{2}$},](v1)(v2)
    \Edge[lw=0.025in,style={color=cv1v3,double=cfv1v3},labelstyle={sloped,pos=0.5,text=clv1v3,},label=\hbox{$x y^{3}$},](v1)(v3)
    \Edge[lw=0.025in,style={color=cv1v4,double=cfv1v4},labelstyle={sloped,pos=0.5,text=clv1v4,},label=\hbox{$x y^{4}$},](v1)(v4)
    \Edge[lw=0.025in,style={color=cv2v3,double=cfv2v3},labelstyle={sloped,left,text=clv2v3,},label=\hbox{$x^{2} y^{3}$},](v2)(v3)
    \Edge[lw=0.025in,style={color=cv2v4,double=cfv2v4},labelstyle={pos=0.5,text=clv2v4,},label=\hbox{$x^{2} y^{4}$},](v2)(v4)
    \Edge[lw=0.1in,style={color=cv3v4,double=cfv3v4},labelstyle={sloped,pos=0.5,text=clv3v4,},label=\hbox{$x^{3} y^{4}$},](v3)(v4)
    \Loop[dist=4.0in,dir=EA,style={color=cv4v4,double=cfv4v4},labelstyle={sloped,below,text=clv4v4,},label=\hbox{$w^{5}$},](v4)
    %
    \end{tikzpicture}

GraphLatex class and functions
------------------------------
"""
# ****************************************************************************
#       Copyright (C) 2009 Robert Beezer <beezer@ups.edu>
#       Copyright (C) 2009 Fidel Barrera Cruz <fidel.barrera@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_function
from sage.misc.latex import latex


def check_tkz_graph():
    r"""
    Check if the proper LaTeX packages for the ``tikzpicture`` environment are
    installed in the user's environment, and issue a warning otherwise.

    The warning is only issued on the first call to this function. So any
    doctest that illustrates the use of the tkz-graph packages should call this
    once as having random output to exhaust the warnings before testing output.

    See also :meth:`sage.misc.latex.Latex.check_file`

    TESTS::

        sage: from sage.graphs.graph_latex import check_tkz_graph
        sage: check_tkz_graph()  # random - depends on TeX installation
        sage: check_tkz_graph()  # at least the second time, so no output
    """
    latex.check_file("tikz.sty", """This package is required to render graphs in LaTeX.
Visit '...'.
""")
    latex.check_file("tkz-graph.sty", """This package is required to render graphs in LaTeX.
Visit 'https://www.ctan.org/pkg/tkz-graph'.
""")
    latex.check_file("tkz-berge.sty", """This package is required to render graphs in LaTeX.
Visit 'https://www.ctan.org/pkg/tkz-berge'.
""")


def have_tkz_graph():
    r"""
    Return ``True`` if the proper LaTeX packages for the ``tikzpicture``
    environment are installed in the user's environment, namely ``tikz``,
    ``tkz-graph`` and ``tkz-berge``.

    The result is cached.

    See also :meth:`sage.misc.latex.Latex.has_file`

    TESTS::

        sage: from sage.graphs.graph_latex import have_tkz_graph
        sage: have_tkz_graph()  # random - depends on TeX installation
        sage: have_tkz_graph() in [True, False]
        True
    """
    return latex.has_file("tikz.sty") and latex.has_file("tkz-graph.sty") and latex.has_file("tkz-berge.sty")


@cached_function
def setup_latex_preamble():
    r"""
    Add appropriate ``\usepackage{...}``, and other instructions to the latex
    preamble for the packages that are needed for processing graphs(``tikz``,
    ``tkz-graph``, ``tkz-berge``), if available in the ``LaTeX`` installation.

    See also :meth:`sage.misc.latex.Latex.add_package_to_preamble_if_available`.

    EXAMPLES::

        sage: sage.graphs.graph_latex.setup_latex_preamble()

    TESTS::

        sage: ("\\usepackage{tikz}" in latex.extra_preamble()) == latex.has_file("tikz.sty")
        True
    """
    latex.add_package_to_preamble_if_available("tikz")
    latex.add_package_to_preamble_if_available("tkz-graph")
    latex.add_package_to_preamble_if_available("tkz-berge")
    if have_tkz_graph():
        latex.add_to_preamble("\\usetikzlibrary{arrows,shapes}")


class GraphLatex(SageObject):
    r"""
    A class to hold, manipulate and employ options for converting
    a graph to LaTeX.

    This class serves two purposes.  First it holds the values of various
    options designed to work with the ``tkz-graph`` LaTeX package for rendering
    graphs.  As such, a graph that uses this class will hold a reference to it.
    Second, this class contains the code to convert a graph into the
    corresponding LaTeX constructs, returning a string.

    EXAMPLES::

        sage: from sage.graphs.graph_latex import GraphLatex
        sage: opts = GraphLatex(graphs.PetersenGraph())
        sage: opts
        LaTeX options for Petersen graph: {}
        sage: g = graphs.PetersenGraph()
        sage: opts = g.latex_options()
        sage: g == loads(dumps(g))
        True
    """

    #  These are the "allowed" options for a graph, private to the class,
    #  along with their default value and description
    #  This allows intelligent errors when non-existent options are referenced
    #  Additionally, for each new option added here:
    #    1.  Document values in GraphLatex.set_option() docstring
    #    2.  Describe also in docstring for the sage.graphs.graph_latex module
    #
    # TODO: use some standard option handling mechanism
    # This dictionary could also contain type information (list of admissible
    # values) and a description
    # See e.g. @option
    __graphlatex_options = {'tkz_style': 'Custom',
                            'format': 'tkz_graph',
                            'layout': 'acyclic',
                            'prog': 'dot',
                            'units': 'cm',
                            'scale': 1.0,
                            'graphic_size': (5, 5),
                            'margins': (0, 0, 0, 0),
                            'vertex_color': 'black',
                            'vertex_colors': {},
                            'vertex_fill_color': 'white',
                            'vertex_fill_colors': {},
                            'vertex_shape': 'circle',
                            'vertex_shapes': {},
                            'vertex_size': 1.0,
                            'vertex_sizes': {},
                            'vertex_labels': True,
                            'vertex_labels_math': True,
                            'vertex_label_color': 'black',
                            'vertex_label_colors': {},
                            'vertex_label_placement': 'center',
                            'vertex_label_placements': {},
                            'edge_options': (),
                            'edge_color': 'black',
                            'edge_colors': {},
                            'edge_fills': False,
                            'edge_fill_color': 'black',
                            'edge_fill_colors': {},
                            'edge_thickness': 0.1,
                            'edge_thicknesses': {},
                            'edge_labels': False,
                            'edge_labels_math': True,
                            'edge_label_color': 'black',
                            'edge_label_colors': {},
                            'edge_label_sloped': True,
                            'edge_label_slopes': {},
                            'edge_label_placement': 0.50,
                            'edge_label_placements': {},
                            'loop_placement': (3.0, 'NO'),
                            'loop_placements': {},
                            'color_by_label': False,
                            'rankdir': 'down',
                            'subgraph_clusters': []}

    def __init__(self, graph, **options):
        r"""
        Return a GraphLatex object, which holds all the parameters needed for
        creating a LaTeX string that will be rendered as a picture of the graph.

        See :mod:`sage.graphs.graph_latex` for more documentation.

        EXAMPLES::

            sage: from sage.graphs.graph_latex import GraphLatex
            sage: GraphLatex(graphs.PetersenGraph())
            LaTeX options for Petersen graph: {}
        """
        self._graph = graph
        self._options = {}
        self.set_options(**options)

    def __eq__(self, other):
        r"""
        Two :class:`sage.graphs.graph_latex.GraphLatex` objects are equal if
        their options are equal.

        The graphs they are associated with are ignored in the comparison.

        TESTS::

            sage: from sage.graphs.graph_latex import GraphLatex
            sage: opts1 = GraphLatex(graphs.PetersenGraph())
            sage: opts2 = GraphLatex(graphs.CompleteGraph(10))
            sage: opts1.set_option('tkz_style', 'Art')
            sage: opts2.set_option('tkz_style', 'Art')
            sage: opts1 == opts2
            True
            sage: opts2.set_option('tkz_style', 'Normal')
            sage: opts1 == opts2
            False
        """
        if not isinstance(other, GraphLatex):
            return False
        return self._options == other._options

    def _repr_(self):
        r"""
        Return a string representation of a
        :class:`sage.graphs.graph_latex.GraphLatex` object which includes the
        name of the graph and the dictionary of current options.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: opts = g.latex_options()
            sage: opts.set_option('tkz_style', 'Classic')
            sage: opts.set_option('vertex_size', 3.6)
            sage: print(opts._repr_())
            LaTeX options for Petersen graph: {'tkz_style': 'Classic', 'vertex_size': 3.60000000000000}
        """
        return "LaTeX options for %s: %s" % (self._graph, self._options)

    def set_option(self, option_name, option_value=None):
        r"""
        Set, modify, clear a LaTeX option for controlling the rendering of a
        graph.

        The possible options are documented here, because ultimately it is this
        routine that sets the values. However, the
        :meth:`sage.graphs.generic_graph.GenericGraph.set_latex_options` method
        is the easiest way to set options, and allows several to be set at once.

        INPUT:

        - ``option_name`` -- string for a latex option contained in the list
          ``sage.graphs.graph_latex.GraphLatex.__graphlatex_options``;
          a :exc:`ValueError` is raised if the option is not allowed

        - ``option_value`` -- a value for the option; if omitted, or set to
          ``None``, the option will use the default value

        The output can be either handled internally by ``Sage``, or delegated to
        the external software ``dot2tex`` and ``graphviz``. This is controlled
        by the option ``format``:

        - ``format`` -- string (default: ``'tkz_graph'``); either ``'dot2tex'``
          or ``'tkz_graph'``

        If format is ``'dot2tex'``, then all the LaTeX generation will be
        delegated to ``dot2tex`` (which must be installed).

        For ``tkz_graph``, the possible option names, and associated values are
        given below. This first group allows you to set a style for a graph and
        specify some sizes related to the eventual image. (For more information
        consult the documentation for the ``tkz-graph`` package.)

        - ``tkz_style`` -- string (default: ``'Custom'``); the name of a
          pre-defined ``tkz-graph`` style such as ``'Shade'``, ``'Art'``,
          ``'Normal'``, ``'Dijkstra'``, ``'Welsh'``, ``'Classic'``, and
          ``'Simple'``, or the string ``'Custom'``.  Using one of these styles
          alone will often give a reasonably good drawing with minimal
          effort. For a custom appearance set this to ``'Custom'`` and use the
          options described below to override the default values.

        - ``units`` -- string (default: ``'cm'``); a natural unit of
          measurement used for all dimensions.  Possible values are: ``'in'``,
          ``'mm'``, ``'cm'``, ``'pt'``, ``'em'``, ``'ex'``.

        - ``scale`` -- float (default: `1.0`); a dimensionless number that
          multiplies every linear dimension. So you can design at sizes you are
          accustomed to, then shrink or expand to meet other needs. Though fonts
          do not scale.

        - ``graphic_size`` -- tuple (default: ``(5, 5)``); overall dimensions
          (width, length) of the bounding box around the entire graphic image

        - ``margins`` -- 4-tuple (default: ``(0, 0, 0, 0)``); portion of graphic
          given over to a plain border as a tuple of four numbers: (left, right,
          top, bottom). These are subtracted from the ``graphic_size`` to
          create the area left for the vertices of the graph itself.  Note that
          the processing done by Sage will trim the graphic down to the minimum
          possible size, removing any border. So this is only useful if you use
          the latex string in a latex document.

        If not using a pre-built style the following options are used, so the
        following defaults will apply.  It is not possible to begin with a
        pre-built style and modify it (other than editing the latex string by
        hand after the fact).

        - ``vertex_color`` -- (default: ``'black'``) a single color to use as
          the default for outline of vertices. For the ``sphere`` shape this
          color is used for the entire vertex, which is drawn with a 3D shading.
          Colors must be specified as a string recognized by the matplotlib
          library: a standard color name like ``'red'``, or a hex string like
          ``'#2D87A7'``, or a single character from the choices ``'rgbcmykw'``.
          Additionally, a number between 0 and 1 will create a grayscale value.
          These color specifications are consistent throughout the options for
          a ``tikzpicture``.

        - ``vertex_colors`` -- dictionary whose keys are vertices of the graph
          and whose values are colors. These will be used to color the outline
          of vertices. See the explanation above for the ``vertex_color`` option
          to see possible values. These values need only be specified for a
          proper subset of the vertices. Specified values will supersede a
          default value.

        - ``vertex_fill_color`` -- (default: ``'white'``) a single color to use
          as the default for the fill color of vertices. See the explanation
          above for the ``vertex_color`` option to see possible values. This
          color is ignored for the ``sphere`` vertex shape.

        - ``vertex_fill_colors`` -- dictionary whose keys are vertices of the
          graph and whose values are colors. These will be used to fill the
          interior of vertices. See the explanation above for the
          ``vertex_color`` option to see possible values. These values need only
          be specified for a proper subset of the vertices. Specified values
          will supersede a default value.

        - ``vertex_shape`` -- string (default: ``'circle'``); specifies the
          shape of the vertices. Allowable values are ``'circle'``,
          ``'sphere'``, ``'rectangle'``, ``'diamond'``. The sphere shape has a
          3D look to its coloring and is uses only one color, that specified by
          ``vertex_color`` and ``vertex_colors``, which are normally used for
          the outline of the vertex.

        - ``vertex_shapes`` -- dictionary whose keys are vertices of the graph
          and whose values are shapes. See ``vertex_shape`` for the allowable
          possibilities.

        - ``vertex_size`` -- float (default: 1.0); the minimum size of a vertex
          as a number. Vertices will expand to contain their labels if the
          labels are placed inside the vertices. If you set this value to zero
          the vertex will be as small as possible (up to tkz-graph's "inner sep"
          parameter), while still containing labels. However, if labels are not
          of a uniform size, then the vertices will not be either.

        - ``vertex_sizes`` -- dictionary of sizes for some of the vertices

        - ``vertex_labels`` -- boolean (default: ``True``); determine whether or
          not to display the vertex labels.  If ``False`` subsequent options
          about vertex labels are ignored.

        - ``vertex_labels_math`` -- boolean (default: ``True``); when ``True``,
          if a label is a string that begins and ends with dollar signs, then
          the string will be rendered as a latex string.  Otherwise, the label
          will be automatically subjected to the ``latex()`` method and rendered
          accordingly. If ``False`` the label is rendered as its textual
          representation according to the ``_repr`` method. Support for
          arbitrarily-complicated mathematics is not especially robust.

        - ``vertex_label_color`` -- (default: ``'black'``) a single color to
          use as the default for labels of vertices. See the explanation above
          for the ``vertex_color`` option to see possible values.

        - ``vertex_label_colors`` -- dictionary whose keys are vertices of the
          graph and whose values are colors. These will be used for the text of
          the labels of vertices. See the explanation above for the
          ``vertex_color`` option to see possible values. These values need only
          be specified for a proper subset of the vertices. Specified values
          will supersede a default value.

        - ``vertex_label_placement`` -- (default: ``'center'``) if ``'center'``
          the label is centered in the interior of the vertex and the vertex
          will expand to contain the label. Giving instead a pair of numbers
          will place the label exterior to the vertex at a certain distance from
          the edge, and at an angle to the positive x-axis, similar in spirit to
          polar coordinates.

        - ``vertex_label_placements`` -- dictionary of placements indexed by
          the vertices. See the explanation for ``vertex_label_placement`` for
          the possible values.

        - ``edge_color`` -- (default: ``'black'``) a single color to use as the
          default for an edge. See the explanation above for the
          ``vertex_color`` option to see possible values.

        - ``edge_colors`` -- dictionary whose keys are edges of the graph and
          whose values are colors. These will be used to color the edges. See
          the explanation above for the ``vertex_color`` option to see possible
          values. These values need only be specified for a proper subset of the
          vertices. Specified values will supersede a default value.

        - ``edge_fills`` -- boolean (default: ``False``); whether an edge has a
          second color running down the middle. This can be a useful effect for
          highlighting edge crossings.

        - ``edge_fill_color`` -- (default: ``'black'``) a single color to use
          as the default for the fill color of an edge. The boolean switch
          ``edge_fills`` must be set to ``True`` for this to have an effect. See
          the explanation above for the ``vertex_color`` option to see possible
          values.

        - ``edge_fill_colors`` -- dictionary whose keys are edges of the graph
          and whose values are colors. See the explanation above for the
          ``vertex_color`` option to see possible values. These values need
          only be specified for a proper subset of the vertices. Specified
          values will supersede a default value.

        - ``edge_thickness`` -- float (default: 0.1); specifies the width of the
          edges. Note that ``tkz-graph`` does not interpret this number for
          loops.

        - ``edge_thicknesses`` -- dictionary of thicknesses for some of the
          edges of a graph. These values need only be specified for a proper
          subset of the vertices. Specified values will supersede a default
          value.

        - ``edge_labels`` -- boolean (default: ``False``); determine if edge
          labels are shown. If ``False`` subsequent options about edge labels
          are ignored.

        - ``edge_labels_math`` -- boolean (default: ``True``); control how edge
          labels are rendered. Read the explanation for the
          ``vertex_labels_math`` option, which behaves identically.  Support for
          arbitrarily-complicated mathematics is not especially robust.

        - ``edge_label_color`` -- (default: ``'black'``) a single color to use
          as the default for labels of edges. See the explanation above for the
          ``vertex_color`` option to see possible values.

        - ``edge_label_colors`` -- dictionary whose keys are edges of the
          graph and whose values are colors. These will be used for the text of
          the labels of edges. See the explanation above for the
          ``vertex_color`` option to see possible values. These values need only
          be specified for a proper subset of the vertices. Specified values
          will supersede a default value. Note that labels must be used for this
          to have any effect, and no care is taken to ensure that label and fill
          colors work well together.

        - ``edge_label_sloped`` -- boolean (default: ``True``); specifies how
          edge labels are place. ``False`` results in a horizontal label, while
          ``True`` means the label is rotated to follow the direction of the
          edge it labels.

        - ``edge_label_slopes`` -- dictionary of booleans, indexed by some
          subset of the edges.  See the ``edge_label_sloped`` option for a
          description of sloped edge labels.

        - ``edge_label_placement`` -- (default: 0.50) either a number between
          0.0 and 1.0, or one of: ``'above'``, ``'below'``, ``'left'``,
          ``'right'``. These adjust the location of an edge label along an
          edge. A number specifies how far along the edge the label is located.
          ``'left'`` and ``'right'`` are conveniences. ``'above'`` and
          ``'below'`` move the label off the edge itself while leaving it near
          the midpoint of the edge. The default value of ``0.50`` places the
          label on the midpoint of the edge.

        - ``edge_label_placements`` -- dictionary of edge placements, indexed
          by the edges.  See the ``edge_label_placement`` option for a
          description of the allowable values.

        - ``loop_placement`` -- (default: ``(3.0, 'NO')``); determine how loops
          are rendered.  the first element of the pair is a distance, which
          determines how big the loop is and the second element is a string
          specifying a compass point (North, South, East, West) as one of
          ``'NO'``, ``'SO'``, ``'EA'``, ``'WE'``.

        - ``loop_placements`` -- dictionary of loop placements.  See the
          ``loop_placements`` option for the allowable values.  While loops are
          technically edges, this dictionary is indexed by vertices.

        For the ``'dot2tex'`` format, the possible option names and associated
        values are given below:

        - ``prog`` -- string; the program used for the layout. It must be a
          string corresponding to one of the software of the graphviz suite:
          ``'dot'``, ``'neato'``, ``'twopi'``, ``'circo'`` or ``'fdp'``.

        - ``edge_labels`` -- boolean (default: ``False)``; whether to display
          the labels on edges

        - ``edge_colors`` -- a color; can be used to set a global color to the
          edge of the graph

        - ``color_by_label`` -- boolean (default: ``False``); colors the edges
          according to their labels

        - ``subgraph_clusters`` -- (default: ``[]``) a list of lists of
          vertices, if supported by the layout engine, nodes belonging to the
          same cluster subgraph are drawn together, with the entire drawing of
          the cluster contained within a bounding rectangle.

        OUTPUT: none; success happens silently

        EXAMPLES:

        Set, then modify, then clear the ``tkz_style`` option, and finally show
        an error for an unrecognized option name::

            sage: g = graphs.PetersenGraph()
            sage: opts = g.latex_options()
            sage: opts
            LaTeX options for Petersen graph: {}
            sage: opts.set_option('tkz_style', 'Art')
            sage: opts
            LaTeX options for Petersen graph: {'tkz_style': 'Art'}
            sage: opts.set_option('tkz_style', 'Simple')
            sage: opts
            LaTeX options for Petersen graph: {'tkz_style': 'Simple'}
            sage: opts.set_option('tkz_style')
            sage: opts
            LaTeX options for Petersen graph: {}
            sage: opts.set_option('bad_name', 'nonsense')
            Traceback (most recent call last):
            ...
            ValueError: bad_name is not a LaTeX option for a graph.

        See :meth:`sage.graphs.generic_graph.GenericGraph.layout_graphviz` for
        installation instructions for ``graphviz`` and ``dot2tex``. Furthermore,
        pgf >= 2.00 should be available inside LaTeX's tree for LaTeX
        compilation (e.g. when using ``view``). In case your LaTeX distribution
        does not provide it, here are short instructions:

           - download pgf from http://sourceforge.net/projects/pgf/
           - unpack it in ``/usr/share/texmf/tex/generic`` (depends on your system)
           - clean out remaining pgf files from older version
           - run texhash

        TESTS:

        These test all of the options and one example of each allowable proper
        input. They should all execute silently. ::

            sage: G = Graph()
            sage: G.add_edge((0,1))
            sage: opts = G.latex_options()
            sage: opts.set_option('tkz_style', 'Custom')
            sage: opts.set_option('tkz_style', 'Art')
            sage: opts.set_option('format', 'tkz_graph')
            sage: opts.set_option('layout', 'acyclic')
            sage: opts.set_option('prog', 'dot')
            sage: opts.set_option('units', 'cm')
            sage: opts.set_option('scale', 1.0)
            sage: opts.set_option('graphic_size', (5, 5))
            sage: opts.set_option('margins', (0,0,0,0))
            sage: opts.set_option('vertex_color', 'black')
            sage: opts.set_option('vertex_colors', {0:'#ABCDEF'})
            sage: opts.set_option('vertex_fill_color', 'white')
            sage: opts.set_option('vertex_fill_colors', {0:'c'})
            sage: opts.set_option('vertex_shape', 'circle')
            sage: opts.set_option('vertex_shapes', {0:'sphere'})
            sage: opts.set_option('vertex_size', 1.0)
            sage: opts.set_option('vertex_sizes', {0:3.4})
            sage: opts.set_option('vertex_labels', True)
            sage: opts.set_option('vertex_labels_math', True)
            sage: opts.set_option('vertex_label_color', 'black')
            sage: opts.set_option('vertex_label_colors', {0:'.23'})
            sage: opts.set_option('vertex_label_placement', 'center')
            sage: opts.set_option('vertex_label_placement', (3, 4.2))
            sage: opts.set_option('vertex_label_placements', {0:'center'})
            sage: opts.set_option('vertex_label_placements', {0:(4.7,1)})
            sage: opts.set_option('edge_color', 'black')
            sage: opts.set_option('edge_colors', {(0,1):'w'})
            sage: opts.set_option('edge_fills', False)
            sage: opts.set_option('edge_fill_color', 'black')
            sage: opts.set_option('edge_fill_colors', {(0,1):"#123456"})
            sage: opts.set_option('edge_thickness', 0.1)
            sage: opts.set_option('edge_thicknesses', {(0,1):5.2})
            sage: opts.set_option('edge_labels', False)
            sage: opts.set_option('edge_labels_math', True)
            sage: opts.set_option('edge_label_color', 'black')
            sage: opts.set_option('edge_label_colors', {(0,1):'red'})
            sage: opts.set_option('edge_label_sloped', True)
            sage: opts.set_option('edge_label_slopes', {(0,1): False})
            sage: opts.set_option('edge_label_placement', 'left')
            sage: opts.set_option('edge_label_placement', 0.50)
            sage: opts.set_option('edge_label_placements', {(0,1):'above'})
            sage: opts.set_option('edge_label_placements', {(0,1):0.75})
            sage: opts.set_option('loop_placement', (3.0, 'NO'))
            sage: opts.set_option('loop_placements', {0:(5.7,'WE')})
            sage: opts.set_option('subgraph_clusters', [[0,1]])

        These test some of the logic of possible failures. Some tests, such as
        inputs of colors, are handled by somewhat general sections of code and
        are not tested for each possible option. ::

            sage: G=Graph()
            sage: G.add_edge((0,1))
            sage: opts = G.latex_options()
            sage: opts.set_option('tkz_style', 'Crazed')
            Traceback (most recent call last):
            ...
            ValueError: tkz_style is not "Custom", nor an implemented tkz-graph style
            sage: opts.set_option('format', 'NonExistent')
            Traceback (most recent call last):
            ...
            ValueError: format option must be one of: tkz_graph, dot2tex not NonExistent
            sage: opts.set_option('units', 'furlongs')
            Traceback (most recent call last):
            ...
            ValueError: units option must be one of: in, mm, cm, pt, em, ex, not furlongs
            sage: opts.set_option('graphic_size', (1,2,3))
            Traceback (most recent call last):
            ...
            ValueError: graphic_size option must be an ordered pair, not (1, 2, 3)
            sage: opts.set_option('margins', (1,2,3))
            Traceback (most recent call last):
            ...
            ValueError: margins option must be 4-tuple, not (1, 2, 3)
            sage: opts.set_option('vertex_color', 'chartruse')
            Traceback (most recent call last):
            ...
            ValueError: vertex_color option needs to be a matplotlib color (always as a string), not chartruse
            sage: opts.set_option('vertex_labels_math', 'maybe')
            Traceback (most recent call last):
            ...
            ValueError: vertex_labels_math option must be True or False, not maybe
            sage: opts.set_option('vertex_shape', 'decagon')
            Traceback (most recent call last):
            ...
            ValueError: vertex_shape option must be the shape of a vertex, not decagon
            sage: opts.set_option('scale', 'big')
            Traceback (most recent call last):
            ...
            ValueError: scale option must be a positive number, not big
            sage: opts.set_option('scale', -6)
            Traceback (most recent call last):
            ...
            ValueError: scale option must be a positive number, not -6
            sage: opts.set_option('vertex_label_placement', (2,-4))
            Traceback (most recent call last):
            ...
            ValueError: vertex_label_placement option must be None, or a pair of positive numbers, not (2, -4)
            sage: opts.set_option('edge_label_placement', 3.6)
            Traceback (most recent call last):
            ...
            ValueError: edge_label_placement option must be a number between 0.0 and 1.0 or a place (like "above"), not 3.60000000000000
            sage: opts.set_option('loop_placement', (5,'SW'))
            Traceback (most recent call last):
            ...
            ValueError: loop_placement option must be a pair that is a positive number followed by a compass point abbreviation, not (5, 'SW')
            sage: opts.set_option('vertex_fill_colors', {0:'#GG0000'})
            Traceback (most recent call last):
            ...
            ValueError: vertex_fill_colors option for 0 needs to be a matplotlib color (always as a string), not #GG0000
            sage: opts.set_option('vertex_sizes', {0:-10})
            Traceback (most recent call last):
            ...
            ValueError: vertex_sizes option for 0 needs to be a positive number, not -10
            sage: opts.set_option('edge_label_slopes', {(0,1):'possibly'})
            Traceback (most recent call last):
            ...
            ValueError: edge_label_slopes option for (0, 1) needs to be True or False, not possibly
            sage: opts.set_option('vertex_shapes', {0:'pentagon'})
            Traceback (most recent call last):
            ...
            ValueError: vertex_shapes option for 0 needs to be a vertex shape, not pentagon
            sage: opts.set_option('vertex_label_placements', {0:(1,2,3)})
            Traceback (most recent call last):
            ...
            ValueError: vertex_label_placements option for 0 needs to be None or a pair of positive numbers, not (1, 2, 3)
            sage: opts.set_option('edge_label_placements', {(0,1):'partway'})
            Traceback (most recent call last):
            ...
            ValueError: edge_label_placements option for (0, 1) needs to be a number between 0.0 and 1.0 or a place (like "above"), not partway
            sage: opts.set_option('loop_placements', {0:(-3,'WE')})
            Traceback (most recent call last):
            ...
            ValueError: loop_placements option for 0 needs to be a positive number and a compass point (like "EA"), not (-3, 'WE')
            sage: opts.set_option('margins', (1,2,3,-5))
            Traceback (most recent call last):
            ...
            ValueError: margins option of (1, 2, 3, -5) cannot contain -5
        """
        # TODO: Needed improvements, possible extensions, dubious ideas
        # - digraph edges should be optionally curved or straight with perhaps a
        #   variable curvature (exit angle from vertex).  Always curved now to
        #   allow for bidirectional.
        # - the "draw" option will make boxes around labels as extensions of the
        #   edge color and thickness
        # - edge labels can have colored backgrounds (which look like fills when
        #   boxed.
        # - edge label fonts can be sized (latex style), which will make scaling
        #   work totally
        # - edges can be dotted or dashed, Beezer suggests calling this "edge
        #   shape" to mirror vertex shapes
        # - "line width" works for vertices, should be configurable
        # - allow injection of latex code to style a pre-built style for
        #   example, \SetUpVertex[style={fill=green}] could override color
        #   selection in a style like "Art"
        # - "inner sep" is distance from vertex label to edge of vertex this
        #   should be set as small as possible - but bigger than the line width.
        # - aspect ratio could be preserved, see hints near creation of affine
        #   transformation.
        # - "outer sep" causes edges to stop some distance before reaching
        #   vertices.  Seems of limited value.
        # - Multi-edges are not supported.  Need to recognize them, twiddle keys
        #   in dictionaries, plot with a spectrum of bends.
        # Seems like a substantial project.

        from matplotlib.colors import ColorConverter
        from sage.rings.integer import Integer
        from sage.rings.real_mpfr import RealLiteral

        cc = ColorConverter()  # used as a color tester

        if option_name not in GraphLatex.__graphlatex_options:
            raise ValueError("%s is not a LaTeX option for a graph." % option_name)
        if option_value is None:    # clear the option, if set
            if option_name in self._options:
                del self._options[option_name]
        else:
            # Test options here when attempt to set
            name = option_name
            value = option_value
            #
            # Tuples of constants
            #
            formats = ('tkz_graph', 'dot2tex')
            styles = ('Custom', 'Shade', 'Art', 'Normal', 'Dijkstra', 'Welsh', 'Classic', 'Simple')
            unit_names = ('in', 'mm', 'cm', 'pt', 'em', 'ex')
            shape_names = ('circle', 'sphere', 'rectangle', 'diamond')
            label_places = ('above', 'below', 'right', 'left')
            compass_points = ('NO', 'SO', 'EA', 'WE')
            number_types = (int, Integer, float, RealLiteral)
            #
            # Options with structurally similar tests
            #
            boolean_options = ('vertex_labels', 'vertex_labels_math', 'edge_fills',
                               'edge_labels', 'edge_labels_math', 'edge_label_sloped')
            color_options = ('vertex_color', 'vertex_fill_color', 'vertex_label_color',
                             'edge_color', 'edge_fill_color', 'edge_label_color')
            color_dicts = ('vertex_colors', 'vertex_fill_colors', 'vertex_label_colors',
                           'edge_colors', 'edge_fill_colors', 'edge_label_colors')
            boolean_dicts = ('edge_label_slopes',)
            positive_scalars = ('scale', 'vertex_size', 'edge_thickness')
            positive_scalar_dicts = ('vertex_sizes', 'edge_thicknesses')
            positive_tuples = ('graphic_size', 'margins')
            #
            #  Checks/test on single values (ie graph-wide defaults)
            #
            if name == 'tkz_style' and value not in styles:
                raise ValueError('%s is not "Custom", nor an implemented tkz-graph style' % name)
            elif name == 'format' and value not in formats:
                raise ValueError('%s option must be one of: tkz_graph, dot2tex not %s' % (name, value))
            elif name == 'units' and value not in unit_names:
                raise ValueError('%s option must be one of: in, mm, cm, pt, em, ex, not %s' % (name, value))
            elif name == 'graphic_size' and not (isinstance(value, tuple) and (len(value) == 2)):
                raise ValueError('%s option must be an ordered pair, not %s' % (name, value))
            elif name == 'margins' and not ((isinstance(value, tuple)) and (len(value) == 4)):
                raise ValueError('%s option must be 4-tuple, not %s' % (name, value))
            elif name in color_options:
                try:
                    cc.to_rgb(value)
                except Exception:
                    raise ValueError('%s option needs to be a matplotlib color (always as a string), not %s' % (name, value))
            elif name in boolean_options and not isinstance(value, bool):
                raise ValueError('%s option must be True or False, not %s' % (name, value))
            elif name == 'vertex_shape' and value not in shape_names:
                raise ValueError('%s option must be the shape of a vertex, not %s' % (name, value))
            elif name in positive_scalars and not (type(value) in number_types and (value >= 0.0)):
                raise ValueError('%s option must be a positive number, not %s' % (name, value))
            elif (name == 'vertex_label_placement' and value != 'center' and
                  not (isinstance(value, tuple) and len(value) == 2 and
                       type(value[0]) in number_types and value[0] >= 0 and
                       type(value[1]) in number_types and value[1] >= 0)):
                raise ValueError('%s option must be None, or a pair of positive numbers, not %s' % (name, value))
            elif (name == 'edge_label_placement' and
                  not ((type(value) in number_types and 0 <= value <= 1)
                       or value in label_places)):
                raise ValueError('%s option must be a number between 0.0 and 1.0 or a place (like "above"), not %s' % (name, value))
            elif (name == 'loop_placement' and
                  not (isinstance(value, tuple) and len(value) == 2 and
                       value[0] >= 0 and value[1] in compass_points)):
                raise ValueError('%s option must be a pair that is a positive number followed by a compass point abbreviation, not %s' % (name, value))
            #
            #  Checks/test on dictionaries of values (ie per-vertex or per-edge defaults)
            #
            elif name in color_dicts:
                if not isinstance(value, dict):
                    raise TypeError('%s option must be a dictionary, not %s' % (name, value))
                else:
                    for key, c in value.items():
                        try:
                            cc.to_rgb(c)
                        except Exception:
                            raise ValueError('%s option for %s needs to be a matplotlib color (always as a string), not %s' % (name, key, c))
            elif name in positive_scalar_dicts:
                if not isinstance(value, dict):
                    raise TypeError('%s option must be a dictionary, not %s' % (name, value))
                else:
                    for key, x in value.items():
                        if not type(x) in [int, Integer, float, RealLiteral] or not x >= 0.0:
                            raise ValueError('%s option for %s needs to be a positive number, not %s' % (name, key, x))
            elif name in boolean_dicts:
                if not isinstance(value, dict):
                    raise TypeError('%s option must be a dictionary, not %s' % (name, value))
                else:
                    for key, b in value.items():
                        if not isinstance(b, bool):
                            raise ValueError('%s option for %s needs to be True or False, not %s' % (name, key, b))
            elif name == 'vertex_shapes':
                if not isinstance(value, dict):
                    raise TypeError('%s option must be a dictionary, not %s' % (name, value))
                else:
                    for key, s in value.items():
                        if s not in shape_names:
                            raise ValueError('%s option for %s needs to be a vertex shape, not %s' % (name, key, s))
            elif name == 'vertex_label_placements':
                if not isinstance(value, dict):
                    raise TypeError('%s option must be a dictionary, not %s' % (name, value))
                else:
                    for key, p in value.items():
                        if (p != 'center' and
                            not (isinstance(p, tuple) and len(p) == 2 and
                                 type(p[0]) in number_types and p[0] >= 0 and
                                 type(p[1]) in number_types and p[1] >= 0)):
                            raise ValueError('%s option for %s needs to be None or a pair of positive numbers, not %s' % (name, key, p))
            elif name == 'edge_label_placements':
                if not isinstance(value, dict):
                    raise TypeError('%s option must be a dictionary, not %s' % (name, value))
                else:
                    for key, p in value.items():
                        if not (type(p) in [float, RealLiteral] and (0 <= p) and (p <= 1)) and (p not in label_places):
                            raise ValueError('%s option for %s needs to be a number between 0.0 and 1.0 or a place (like "above"), not %s' % (name, key, p))
            elif name == 'loop_placements':
                if not isinstance(value, dict):
                    raise TypeError('%s option must be a dictionary, not %s' % (name, value))
                else:
                    for key, p in value.items():
                        if not ((isinstance(p, tuple)) and (len(p) == 2) and (p[0] >= 0) and (p[1] in compass_points)):
                            raise ValueError('%s option for %s needs to be a positive number and a compass point (like "EA"), not %s' % (name, key, p))
            # These have been verified as tuples before going into this next check
            elif name in positive_tuples:
                for x in value:
                    if type(x) not in [int, Integer, float, RealLiteral] or not x >= 0.0:
                        raise ValueError('%s option of %s cannot contain %s' % (name, value, x))
            #
            # Verified.  Set it.
            self._options[option_name] = option_value

    def set_options(self, **kwds):
        r"""
        Set several LaTeX options for a graph all at once.

        INPUT:

        - ``kwds`` -- any number of option/value pairs to set many graph latex
          options at once (a variable number, in any order). Existing values are
          overwritten, new values are added.  Existing values can be cleared by
          setting the value to ``None``. Errors are raised in the
          :func:`set_option` method.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: opts = g.latex_options()
            sage: opts.set_options(tkz_style='Welsh')
            sage: opts.get_option('tkz_style')
            'Welsh'
        """
        if kwds:
            for name, value in kwds.items():
                self.set_option(name, value)

    def get_option(self, option_name):
        r"""
        Return the current value of the named option.

        INPUT:

        - ``option_name`` -- the name of an option

        OUTPUT:

        If the name is not present in ``__graphlatex_options`` it is an error to
        ask for it. If an option has not been set then the default value is
        returned. Otherwise, the value of the option is returned.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: opts = g.latex_options()
            sage: opts.set_option('tkz_style', 'Art')
            sage: opts.get_option('tkz_style')
            'Art'
            sage: opts.set_option('tkz_style')
            sage: opts.get_option('tkz_style') == "Custom"
            True
            sage: opts.get_option('bad_name')
            Traceback (most recent call last):
            ...
            ValueError: bad_name is not a Latex option for a graph.
        """
        if option_name not in GraphLatex.__graphlatex_options:
            raise ValueError("%s is not a Latex option for a graph." % option_name)
        else:
            if option_name in self._options:
                return self._options[option_name]
            else:
                return GraphLatex.__graphlatex_options[option_name]

    def latex(self):
        r"""
        Return a string in LaTeX representing a graph.

        This is the command that is invoked by
        ``sage.graphs.generic_graph.GenericGraph._latex_`` for a graph, so it
        returns a string of LaTeX commands that can be incorporated into a LaTeX
        document unmodified. The exact contents of this string are influenced by
        the options set via the methods
        :meth:`sage.graphs.generic_graph.GenericGraph.set_latex_options`,
        :meth:`set_option`, and :meth:`set_options`.

        By setting the ``format`` option different packages can be used to
        create the latex version of a graph. Supported packages are
        ``tkz-graph`` and ``dot2tex``.

        EXAMPLES::

            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on TeX installation
            sage: g = graphs.CompleteGraph(2)
            sage: opts = g.latex_options()
            sage: print(opts.latex())
            \begin{tikzpicture}
            \definecolor{cv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv0}{rgb}{1.0,1.0,1.0}
            \definecolor{clv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv1}{rgb}{1.0,1.0,1.0}
            \definecolor{clv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cv0v1}{rgb}{0.0,0.0,0.0}
            %
            \Vertex[style={minimum size=1.0cm,draw=cv0,fill=cfv0,text=clv0,shape=circle},LabelOut=false,L=\hbox{$0$},x=2.5cm,y=5.0cm]{v0}
            \Vertex[style={minimum size=1.0cm,draw=cv1,fill=cfv1,text=clv1,shape=circle},LabelOut=false,L=\hbox{$1$},x=2.5cm,y=0.0cm]{v1}
            %
            \Edge[lw=0.1cm,style={color=cv0v1,},](v0)(v1)
            %
            \end{tikzpicture}

        We check that :issue:`22070` is fixed::

            sage: edges = [(i,(i+1)%3,a) for i,a in enumerate('abc')]
            sage: G_with_labels = DiGraph(edges)
            sage: C = [[0,1], [2]]
            sage: kwds = dict(subgraph_clusters=C,color_by_label=True,prog='dot',format='dot2tex')
            sage: opts = G_with_labels.latex_options()
            sage: opts.set_options(edge_labels=True, **kwds)  # optional - dot2tex graphviz
            sage: latex(G_with_labels)                        # optional - dot2tex graphviz
            \begin{tikzpicture}[>=latex,line join=bevel,]
            %%
            \begin{scope}
              \pgfsetstrokecolor{black}
              \definecolor{strokecol}{rgb}{...};
              \pgfsetstrokecolor{strokecol}
              \definecolor{fillcol}{rgb}{...};
              \pgfsetfillcolor{fillcol}
              \filldraw ... cycle;
            \end{scope}
            \begin{scope}
              \pgfsetstrokecolor{black}
              \definecolor{strokecol}{rgb}{...};
              \pgfsetstrokecolor{strokecol}
              \definecolor{fillcol}{rgb}{...};
              \pgfsetfillcolor{fillcol}
              \filldraw ... cycle;
            \end{scope}
            ...
            \end{tikzpicture}
        """
        format = self.get_option('format')
        if format == "tkz_graph":
            return self.tkz_picture()
        elif format == "dot2tex":
            return self.dot2tex_picture()

    def dot2tex_picture(self):
        r"""
        Call ``dot2tex`` to construct a string of LaTeX commands representing a
        graph as a ``tikzpicture``.

        EXAMPLES::

            sage: g = digraphs.ButterflyGraph(1)
            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on TeX installation
            sage: print(g.latex_options().dot2tex_picture())  # optional - dot2tex graphviz
            \begin{tikzpicture}[>=latex,line join=bevel,]
            %%
              \node (node_...) at (...bp,...bp) [draw,draw=none] {$\left(...\right)$};
              \node (node_...) at (...bp,...bp) [draw,draw=none] {$\left(...\right)$};
              \node (node_...) at (...bp,...bp) [draw,draw=none] {$\left(...\right)$};
              \node (node_...) at (...bp,...bp) [draw,draw=none] {$\left(...\right)$};
              \draw [black,->] (node_...) ..controls (...bp,...bp) and (...bp,...bp)  .. (node_...);
              \draw [black,->] (node_...) ..controls (...bp,...bp) and (...bp,...bp)  .. (node_...);
              \draw [black,->] (node_...) ..controls (...bp,...bp) and (...bp,...bp)  .. (node_...);
              \draw [black,->] (node_...) ..controls (...bp,...bp) and (...bp,...bp)  .. (node_...);
            %
            \end{tikzpicture}

        We make sure :issue:`13624` is fixed::

            sage: G = DiGraph()
            sage: G.add_edge(3333, 88, 'my_label')
            sage: G.set_latex_options(edge_labels=True)
            sage: print(G.latex_options().dot2tex_picture())  # optional - dot2tex graphviz
            \begin{tikzpicture}[>=latex,line join=bevel,]
            %%
            \node (node_...) at (...bp,...bp) [draw,draw=none] {$...$};
              \node (node_...) at (...bp,...bp) [draw,draw=none] {$...$};
              \draw [black,->] (node_...) ..controls (...bp,...bp) and (...bp,...bp)  .. (node_...);
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw (...bp,...bp) node {$\text{\texttt{my{\char`\_}label}}$};
            %
            \end{tikzpicture}

        Check that :issue:`25120` is fixed::

            sage: G = Graph([(0,1)])
            sage: G.set_latex_options(edge_colors = {(0,1): 'red'})
            sage: print(G.latex_options().dot2tex_picture())  # optional - dot2tex graphviz
            \begin{tikzpicture}[>=latex,line join=bevel,]
            ...
            \draw [red,] (node_0) ... (node_1);
            ...
            \end{tikzpicture}

        .. NOTE::

            There is a lot of overlap between what ``tkz_picture`` and
            ``dot2tex`` do. It would be best to merge them! ``dot2tex`` probably
            can work without ``graphviz`` if layout information is provided.
        """
        from sage.graphs.dot2tex_utils import assert_have_dot2tex
        assert_have_dot2tex()

        options = self.__graphlatex_options.copy()
        options.update(self._options)

        # NOTE: the edge_labels option for graphviz_string is opposite
        # see https://github.com/sagemath/sage/issues/25120
        if 'edge_colors' in options:
            edge_colors = options['edge_colors']
            new_edge_colors = {}
            for edge, col in edge_colors.items():
                if col in new_edge_colors:
                    new_edge_colors[col].append(edge)
                else:
                    new_edge_colors[col] = [edge]
            options['edge_colors'] = new_edge_colors

        dotdata = self._graph.graphviz_string(labels='latex', **options)
        import dot2tex
        return dot2tex.dot2tex(dotdata,
                               format='tikz',
                               autosize=True,
                               crop=True,
                               figonly='True',
                               prog=self.get_option('prog')).strip()
        # usepdflatex = True, debug = True)

    def tkz_picture(self):
        r"""
        Return a string of LaTeX commands representing a graph as a
        ``tikzpicture``.

        This routine interprets the graph's properties and the options in
        ``_options`` to render the graph with commands from the ``tkz-graph``
        LaTeX package.

        This requires that the LaTeX optional packages ``tkz-graph`` and
        ``tkz-berge`` be installed. You may also need a current version of the
        pgf package.  If the ``tkz-graph`` and ``tkz-berge`` packages are
        present in the system's TeX installation, the appropriate
        ``\\usepackage{}`` commands will be added to the LaTeX preamble as part
        of the initialization of the graph. If these two packages are not
        present, then this command will return a warning on its first use, but
        will return a string that could be used elsewhere, such as a LaTeX
        document.

        For more information about tkz-graph you can visit
        https://www.ctan.org/pkg/tkz-graph.

        EXAMPLES:

        With a pre-built ``tkz-graph`` style specified, the latex representation
        will be relatively simple. ::

            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on TeX installation
            sage: g = graphs.CompleteGraph(3)
            sage: opts = g.latex_options()
            sage: g.set_latex_options(tkz_style='Art')
            sage: print(opts.tkz_picture())
            \begin{tikzpicture}
            \GraphInit[vstyle=Art]
            %
            \Vertex[L=\hbox{$0$},x=2.5cm,y=5.0cm]{v0}
            \Vertex[L=\hbox{$1$},x=0.0cm,y=0.0cm]{v1}
            \Vertex[L=\hbox{$2$},x=5.0cm,y=0.0cm]{v2}
            %
            \Edge[](v0)(v1)
            \Edge[](v0)(v2)
            \Edge[](v1)(v2)
            %
            \end{tikzpicture}

        Setting the style to "Custom" results in various configurable aspects
        set to the defaults, so the string is more involved. ::

            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on TeX installation
            sage: g = graphs.CompleteGraph(3)
            sage: opts = g.latex_options()
            sage: g.set_latex_options(tkz_style='Custom')
            sage: print(opts.tkz_picture())
            \begin{tikzpicture}
            \definecolor{cv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv0}{rgb}{1.0,1.0,1.0}
            \definecolor{clv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv1}{rgb}{1.0,1.0,1.0}
            \definecolor{clv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cv2}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv2}{rgb}{1.0,1.0,1.0}
            \definecolor{clv2}{rgb}{0.0,0.0,0.0}
            \definecolor{cv0v1}{rgb}{0.0,0.0,0.0}
            \definecolor{cv0v2}{rgb}{0.0,0.0,0.0}
            \definecolor{cv1v2}{rgb}{0.0,0.0,0.0}
            %
            \Vertex[style={minimum size=1.0cm,draw=cv0,fill=cfv0,text=clv0,shape=circle},LabelOut=false,L=\hbox{$0$},x=2.5cm,y=5.0cm]{v0}
            \Vertex[style={minimum size=1.0cm,draw=cv1,fill=cfv1,text=clv1,shape=circle},LabelOut=false,L=\hbox{$1$},x=0.0cm,y=0.0cm]{v1}
            \Vertex[style={minimum size=1.0cm,draw=cv2,fill=cfv2,text=clv2,shape=circle},LabelOut=false,L=\hbox{$2$},x=5.0cm,y=0.0cm]{v2}
            %
            \Edge[lw=0.1cm,style={color=cv0v1,},](v0)(v1)
            \Edge[lw=0.1cm,style={color=cv0v2,},](v0)(v2)
            \Edge[lw=0.1cm,style={color=cv1v2,},](v1)(v2)
            %
            \end{tikzpicture}

        See the introduction to the :mod:`~sage.graphs.graph_latex` module for
        more information on the use of this routine.

        TESTS:

        Graphs with preset layouts that are vertical or horizontal can cause
        problems. First test is a horizontal layout on a path with three
        vertices. ::

            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on TeX installation
            sage: g = graphs.PathGraph(3)
            sage: opts = g.latex_options()
            sage: print(opts.tkz_picture())
            \begin{tikzpicture}
            ...
            \end{tikzpicture}

        Scaling to a bounding box is problematic for graphs with just one
        vertex, or none. ::

            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on TeX installation
            sage: g = graphs.CompleteGraph(1)
            sage: opts = g.latex_options()
            sage: print(opts.tkz_picture())
            \begin{tikzpicture}
            ...
            \end{tikzpicture}

        With the empty graph, an empty tikzfigure is output. ::

            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on TeX installation
            sage: g = Graph()
            sage: opts = g.latex_options()
            sage: print(opts.tkz_picture())
            \begin{tikzpicture}
            %
            %
            %
            \end{tikzpicture}
        """
        # This routine does not handle multiple edges
        # It will properly handle digraphs where a pair of vertices has an edge
        # in each direction, since edges of a digraph are curved.
        if self._graph.has_multiple_edges():
            raise NotImplementedError('it is not possible create a tkz-graph version of a graph with multiple edges')

        from matplotlib.colors import ColorConverter
        from sage.misc.latex import latex

        # On first use of this method, the next call may print warnings as a
        # side effect, but will be silent on any subsequent use.
        check_tkz_graph()

        # Overhead
        cc = ColorConverter()  # .to_rgb method to convert "colors" to triples
        prefix = 'v'  # leading string on internal (to tkz-graph) vertex names

        ####################
        #  Pre-built styles
        ####################

        # We preserve the pre-built style OR get defaults for each option, but
        # we do not mix the two
        style = self.get_option('tkz_style')
        customized = (style == 'Custom')
        # We don't do much for a pre-built style
        # Layout information from the graph
        # And vertex labels (if used) are the latex representation of Sage objects
        if not customized:
            vertex_labels_math = True

        ###################################
        #  Layout, image sizing placement
        ###################################

        units = self.get_option('units')
        scale = self.get_option('scale')
        graphic_size = self.get_option('graphic_size')
        margins = self.get_option('margins')

        # The positions of the vertices will get scaled to fill the specified
        # size of the image, as given by graphic_size.  But first a border is
        # subtracted away and the graph is scaled to fit there.

        # Lower left, upper right corners of box inside borders
        llx = margins[0]
        lly = margins[3]
        urx = graphic_size[0] - margins[1]
        ury = graphic_size[1] - margins[2]
        # width and height of space
        w = urx - llx
        h = ury - lly

        # TODO: Could use self._graph._layout_bounding_box(pos)
        # trans = lambda x,y: [x[0]-y[0],x[1]-y[1]]
        # Determine the spread in the x and y directions (i.e. xmax, ymax)
        # Needs care for perfectly horizontal and vertical layouts

        # We grab the graph's layout (or it is computed as a consequence of the
        # request)
        pos = self._graph.layout()
        if pos.values():
            xmin = min(i[0] for i in pos.values())
            ymin = min(i[1] for i in pos.values())
            xmax = max(i[0] for i in pos.values())
            ymax = max(i[1] for i in pos.values())
        else:
            xmin, xmax, ymin, ymax = 0, 0, 0, 0

        # Linear scaling factors that will be used to scale the image to fit
        # into the bordered region.  Purely horizontal, or purely vertical,
        # layouts get put in the middle of the bounding box by setting the
        # scaling to a constant value on a midline
        xspread = xmax - xmin
        if xspread == 0:
            x_scale = 0.0
            llx = llx + 0.5 * w
        else:
            x_scale = float(w) / xspread
        yspread = ymax - ymin
        if yspread == 0:
            y_scale = 0.0
            lly = lly + 0.5 * h
        else:
            y_scale = float(h) / yspread
        # Could preserve aspect ratio here by setting both scale factors to the
        # minimum and doing a shift of the larger to center
        # A linear function will map layout positions into the bordered graphic
        # space

        def translate(p):
            return ((p[0] - xmin) * x_scale + llx,
                    (p[1] - ymin) * y_scale + lly)

        #############
        #  Vertices
        #############

        # We record the index of each vertex in the graph's list of vertices
        # Which is just a convenience for forming vertex names internal to
        # tkz-graph
        index_of_vertex = {}
        vertex_list = self._graph.vertices(sort=False)
        for u in self._graph:
            index_of_vertex[u] = vertex_list.index(u)

        # Vertex labels can be switched on/off, and we don't record or use this
        # type of extra information if they are switched off
        vertex_labels = self.get_option('vertex_labels')

        # We collect options for vertices, default values and for-some-vertices
        # information. These are combined into dictionaries on a per-vertex
        # basis, for all vertices. This only applies for a custom style
        #
        # Defaults
        #
        if customized:
            dvc = cc.to_rgb(self.get_option('vertex_color'))
            dvfc = cc.to_rgb(self.get_option('vertex_fill_color'))
            dsh = self.get_option('vertex_shape')
            dvs = self.get_option('vertex_size')
            #
            # Default label information, if using vertex labels
            #
            if vertex_labels:
                vertex_labels_math = self.get_option('vertex_labels_math')
                dvlc = cc.to_rgb(self.get_option('vertex_label_color'))
                dvlp = self.get_option('vertex_label_placement')
                # needs test for a pair of numbers, angle and distance (or None)

            # Retrieve dictionaries for selected vertices
            vertex_colors = self.get_option('vertex_colors')
            vertex_fill_colors = self.get_option('vertex_fill_colors')
            vertex_shapes = self.get_option('vertex_shapes')
            vertex_sizes = self.get_option('vertex_sizes')
            if vertex_labels:
                vertex_label_colors = self.get_option('vertex_label_colors')
                vertex_label_placements = self.get_option('vertex_label_placements')

            # Form dictionaries, each indexed for all vertices
            v_color = {}
            vf_color = {}
            v_shape = {}
            v_size = {}
            if vertex_labels:
                vl_color = {}
                vl_placement = {}
            for u in vertex_list:

                c = dvc
                if u in vertex_colors:
                    c = cc.to_rgb(vertex_colors[u])
                v_color[u] = c

                c = dvfc
                if u in vertex_fill_colors:
                    c = cc.to_rgb(vertex_fill_colors[u])
                vf_color[u] = c

                sh = dsh
                if u in vertex_shapes:
                    sh = vertex_shapes[u]
                v_shape[u] = sh

                vs = dvs
                if u in vertex_sizes:
                    vs = vertex_sizes[u]
                v_size[u] = vs

                if vertex_labels:

                    c = dvlc
                    if u in vertex_label_colors:
                        c = cc.to_rgb(vertex_label_colors[u])
                    vl_color[u] = c

                    vlp = dvlp
                    if u in vertex_label_placements:
                        vlp = vertex_label_placements[u]
                        # test vlp here
                    vl_placement[u] = vlp

        is_directed = self._graph.is_directed()

        ##########
        #  Edges
        ##########

        if customized:
            # An "edge fill" is a bit unusual, so we allow it to be turned off
            # as the default.
            edge_fills = self.get_option('edge_fills')

            # Edge labels can be switched on/off, and we don't record or use
            # this type of extra information if they are switched off
            edge_labels = self.get_option('edge_labels')

            # We collect options for edges, default values and for-some-edges
            # information.  These are combined into dictionaries on a per-edge
            # basis, for all edges

            # Defaults

            dec = cc.to_rgb(self.get_option('edge_color'))
            if edge_fills:
                defc = cc.to_rgb(self.get_option('edge_fill_color'))
            det = self.get_option('edge_thickness')

            if edge_labels:
                edge_labels_math = self.get_option('edge_labels_math')
                delc = cc.to_rgb(self.get_option('edge_label_color'))
                dels = self.get_option('edge_label_sloped')
                delp = self.get_option('edge_label_placement')

            # Retrieve dictionaries for selected edges
            edge_colors = self.get_option('edge_colors')
            if edge_fills:
                edge_fill_colors = self.get_option('edge_fill_colors')
            edge_thicknesses = self.get_option('edge_thicknesses')
            if edge_labels:
                edge_label_colors = self.get_option('edge_label_colors')
                edge_label_slopes = self.get_option('edge_label_slopes')
                edge_label_placements = self.get_option('edge_label_placements')

            # Form dictionaries, each indexed for all edges
            #
            # A key of a dictionary indexed by edges may be set for an edge of
            # an undirected graph in the "wrong" order, so we use a "reverse" to
            # test for this case.  Everything formed here conforms to the order
            # used in the graph.

            e_color = {}
            if edge_fills:
                ef_color = {}
            e_thick = {}
            if edge_labels:
                el_color = {}
                el_slope = {}
                el_placement = {}

            for e in self._graph.edges(sort=False):
                edge = (e[0], e[1])
                reverse = (e[1], e[0])

                c = dec
                if edge in edge_colors or (not is_directed and reverse in edge_colors):
                    if edge in edge_colors:
                        c = cc.to_rgb(edge_colors[edge])
                    else:
                        c = cc.to_rgb(edge_colors[reverse])
                e_color[edge] = c

                if edge_fills:
                    c = defc
                    if edge in edge_fill_colors or (not is_directed and reverse in edge_fill_colors):
                        if edge in edge_colors:
                            c = cc.to_rgb(edge_fill_colors[edge])
                        else:
                            c = cc.to_rgb(edge_fill_colors[reverse])
                    ef_color[edge] = c

                et = det
                if edge in edge_thicknesses or (not is_directed and reverse in edge_thicknesses):
                    if edge in edge_thicknesses:
                        et = edge_thicknesses[edge]
                    else:
                        et = edge_thicknesses[reverse]
                e_thick[edge] = et

                if edge_labels:
                    c = delc
                    if edge in edge_label_colors or (not is_directed and reverse in edge_label_colors):
                        if edge in edge_label_colors:
                            c = cc.to_rgb(edge_label_colors[edge])
                        else:
                            c = cc.to_rgb(edge_label_colors[reverse])
                    el_color[edge] = c

                    els = dels
                    if edge in edge_label_slopes or (not is_directed and reverse in edge_label_slopes):
                        if edge in edge_label_slopes:
                            els = edge_label_slopes[edge]
                        else:
                            els = edge_label_slopes[reverse]
                    el_slope[edge] = els

                    elp = delp
                    if edge in edge_label_placements or (not is_directed and reverse in edge_label_placements):
                        if edge in edge_label_placements:
                            elp = edge_label_placements[edge]
                        else:
                            elp = edge_label_placements[reverse]
                    el_placement[edge] = elp

        ##########
        #  Loops
        ##########

        # Loops can be styled much like any other edge by indexing on a pair of
        # two equal vertices though edge thickness is not implemented in
        # tkz-graph!  Size and direction are unique, and are indexed by the
        # vertex rather than on edges.

        # Loop placements are pairs of  length, compass-point
        if customized:
            if self._graph.has_loops():
                dlp = self.get_option('loop_placement')
                loop_placements = self.get_option('loop_placements')
                lp_placement = {}
                for u in vertex_list:
                    lp = dlp
                    if u in loop_placements:
                        lp = loop_placements[u]
                    lp_placement[u] = lp

        ############################
        #  Build the output string
        ############################

        # s is the eventual tkz string
        # Everything should now be in place
        # We build a list and then concatenate it as the return value
        s = ['\\begin{tikzpicture}\n']

        if not customized:
            s += ['\\GraphInit[vstyle=', style, ']\n%\n']

        # Internal strings representing colors are defined here in custom style
        if customized:
            # Define all the colors for the vertices: perimeter, fill, label
            vertex_color_names = {}
            vertex_fill_color_names = {}
            vertex_label_color_names = {}
            for u in vertex_list:
                vertex_color_names[u] = 'c' + prefix + str(index_of_vertex[u])
                s += [r'\definecolor{', vertex_color_names[u], '}{rgb}', '{']
                s += [str(round(v_color[u][0], 4)), ',']
                s += [str(round(v_color[u][1], 4)), ',']
                s += [str(round(v_color[u][2], 4)), '}\n']
                vertex_fill_color_names[u] = 'cf' + prefix + str(index_of_vertex[u])
                s += [r'\definecolor{', vertex_fill_color_names[u], '}{rgb}', '{']
                s += [str(round(vf_color[u][0], 4)), ',']
                s += [str(round(vf_color[u][1], 4)), ',']
                s += [str(round(vf_color[u][2], 4)), '}\n']
                if vertex_labels:
                    vertex_label_color_names[u] = 'cl' + prefix + str(index_of_vertex[u])
                    s += [r'\definecolor{', vertex_label_color_names[u], '}{rgb}{']
                    s += [str(round(vl_color[u][0], 4)), ',']
                    s += [str(round(vl_color[u][1], 4)), ',']
                    s += [str(round(vl_color[u][2], 4)), '}\n']
            # Define all the colors for the edges: perimeter, fill, label
            edge_color_names = {}
            edge_fill_color_names = {}
            edge_label_color_names = {}
            for e in self._graph.edge_iterator():
                edge = (e[0], e[1])
                edge_color_names[edge] = 'c' + prefix + str(index_of_vertex[edge[0]]) + prefix + str(index_of_vertex[edge[1]])
                s += [r'\definecolor{', edge_color_names[edge], '}{rgb}{']
                s += [str(round(e_color[edge][0], 4)), ',']
                s += [str(round(e_color[edge][1], 4)), ',']
                s += [str(round(e_color[edge][2], 4)), '}\n']
                if edge_fills:
                    edge_fill_color_names[edge] = 'cf' + prefix + str(index_of_vertex[edge[0]]) + prefix + str(index_of_vertex[edge[1]])
                    s += [r'\definecolor{', edge_fill_color_names[edge], '}{rgb}{']
                    s += [str(round(ef_color[edge][0], 4)), ',']
                    s += [str(round(ef_color[edge][1], 4)), ',']
                    s += [str(round(ef_color[edge][2], 4)), '}\n']
                if edge_labels:
                    edge_label_color_names[edge] = 'cl' + prefix + str(index_of_vertex[edge[0]]) + prefix + str(index_of_vertex[edge[1]])
                    s += [r'\definecolor{', edge_label_color_names[edge], '}{rgb}{']
                    s += [str(round(el_color[edge][0], 4)), ',']
                    s += [str(round(el_color[edge][1], 4)), ',']
                    s += [str(round(el_color[edge][2], 4)), '}\n']
            s += ['%\n']

        # Create each vertex
        for u in vertex_list:
            s += ['\\Vertex[']
            # colors, shapes, sizes, labels/placement for 'Custom' style
            if customized:
                s += ['style={']  # begin style list
                s += ['minimum size=', str(round(float(scale * v_size[u]), 4)),
                      units, ',']
                s += ['draw=', vertex_color_names[u], ',']
                s += ['fill=', vertex_fill_color_names[u], ',']
                if vertex_labels:
                    s += ['text=', vertex_label_color_names[u], ',']
                if v_shape[u] == 'sphere':
                    s += ['shape=circle,shading=ball,line width=0pt,ball color=', vertex_color_names[u], ',']
                else:
                    s += ['shape=', v_shape[u]]
                s += ['},']  # end style list
                if vertex_labels:
                    if vl_placement[u] == 'center':
                        s += ['LabelOut=false,']
                    else:
                        s += ['LabelOut=true,']
                        s += ['Ldist=', str(round(float(scale * vl_placement[u][0]), 4)), units, ',']
                        s += ['Lpos=', str(round(float(vl_placement[u][1]), 4)), ',']  # degrees, no units
                else:
                    s += ['NoLabel,']
            # vertex label information is available to all pre-built styles
            # but may be ignored by the style, so not apparent
            if vertex_labels or not customized:
                if vertex_labels_math and not (isinstance(u, str) and u[0] == '$' and u[-1] == '$'):
                    lab = r'\hbox{$%s$}' % latex(u)
                else:
                    lab = r'\hbox{%s}' % u
                s += ['L=', lab, ',']
            scaled_pos = translate(pos[u])
            s += ['x=', str(round(float(scale * scaled_pos[0]), 4)), units, ',']
            s += ['y=', str(round(float(scale * scaled_pos[1]), 4)), units]
            s += [']']
            s += ['{', prefix, str(index_of_vertex[u]), '}\n']
        s += ['%\n']

        # Create each edge or loop
        for e in self._graph.edges(sort=False):
            edge = (e[0], e[1])
            loop = e[0] == e[1]
            if loop:
                u = e[0]
                s += ['\\Loop[']
                if customized:
                    s += ['dist=', str(round(float(scale * lp_placement[u][0]), 4)), units, ',']
                    s += ['dir=', lp_placement[u][1], ',']
            else:
                s += ['\\Edge[']
            # colors, shapes, sizes, labels/placement for 'Custom' style
            if customized:
                if not loop:  # lw not available for loops!
                    s += ['lw=', str(round(float(scale * e_thick[edge]), 4)), units, ',']
                s += ['style={']  # begin style list
                if is_directed and not loop:
                    s += ['post, bend right', ',']
                s += ['color=', edge_color_names[edge], ',']
                if edge_fills:
                    s += ['double=', edge_fill_color_names[edge]]
                s += ['},']     # end style list
                if edge_labels:
                    s += ['labelstyle={']
                    if el_slope[edge]:
                        s += ['sloped,']
                    if isinstance(el_placement[edge], str):
                        s += [el_placement[edge], ',']
                    else:
                        s += ['pos=', str(round(float(el_placement[edge]), 4)), ',']  # no units needed
                    s += ['text=', edge_label_color_names[edge], ',']
                    s += ['},']
                    el = self._graph.edge_label(edge[0], edge[1])
                    if edge_labels_math and not (isinstance(el, str) and el[0] == '$' and el[-1] == '$'):
                        lab = r'\hbox{$%s$}' % latex(el)
                    else:
                        lab = r'\hbox{%s}' % el
                    s += ['label=', lab, ',']
            s += [']']
            if not loop:
                s += ['(', prefix, str(index_of_vertex[e[0]]), ')']
            s += ['(', prefix, str(index_of_vertex[e[1]]), ')\n']

        # Wrap it up
        s += ['%\n']
        s += ['\\end{tikzpicture}']

        return ''.join(s)
