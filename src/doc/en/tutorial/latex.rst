*********************************
Sage, LaTeX and Friends
*********************************

AUTHOR:  Rob Beezer (2010-05-23)

Sage and the LaTeX dialect of TeX have an
intensely synergistic relationship. This section aims to
introduce the variety of interactions, beginning with the most
basic and proceeding to the more unusual and arcane.  (So you may
not want to read this entire section on your first pass through
this tutorial.)

Overview
========

It may be easiest to understand the various uses of LaTeX with a
brief overview of the mechanics of the three principal methods
employed by Sage.

    #. Every "object" in Sage is required to have a LaTeX representation.
       You can access this representation by executing, in the notebook or
       at the sage command line, ``latex(foo)`` where ``foo`` is some object
       in Sage.  The output is a string that should render a reasonably accurate
       representation of ``foo`` when used in TeX's math-mode (for example,
       when enclosed between a pair of single dollar signs).  Some examples of
       this follow below.

       In this way, Sage can be used effectively for constructing portions of
       a LaTeX document: create or compute an object in Sage, print ``latex()``
       of the object and cut/paste it into your document.

    #. The notebook interface is configured to use
       `MathJax <http://www.mathjax.org>`_
       to render mathematics
       cleanly in a web browser.  MathJax is an open source JavaScript
       display engine for mathematics that works in all modern
       browsers.  It is able to render a large, but not totally
       complete, subset of TeX.  It has no support for
       things like complicated tables, sectioning or document
       management, as it is oriented towards accurately rendering
       "snippets" of TeX. Seemingly automatic rendering of math in the
       notebook is provided by converting the ``latex()``
       representation of an object (as described above) into a form of
       HTML palatable to MathJax.

       Since MathJax uses its own scalable fonts, it is superior to other methods that
       rely on converting equations, or other snippets of TeX, into static inline images.

    #. At the Sage command-line, or in the notebook when LaTeX code is
       more involved than MathJax can handle, a system-wide installation of
       LaTeX can be employed.  Sage includes almost everything you need to
       build and use Sage, but a significant exception is TeX itself.  So in these
       situations you need to have TeX installed, along with some associated
       conversion utilities, to utilize the full power.

Here we demonstrate some basic uses of the ``latex()`` function. ::

    sage: var('z')
    z
    sage: latex(z^12)
    z^{12}
    sage: latex(integrate(z^4, z))
    \frac{1}{5} \, z^{5}
    sage: latex('a string')
    \text{\texttt{a{ }string}}
    sage: latex(QQ)
    \Bold{Q}
    sage: latex(matrix(QQ, 2, 3, [[2,4,6],[-1,-1,-1]]))
    \left(\begin{array}{rrr}
    2 & 4 & 6 \\
    -1 & -1 & -1
    \end{array}\right)

Basic MathJax functionality is largely automatic in the notebook, but
we can partially demonstrate this support with the ``MathJax`` class.
The ``eval`` function of this class converts a Sage object to its
LaTeX representation and then wraps it in HTML that invokes the CSS
"math" class, which then employs MathJax.  ::

    sage: from sage.misc.latex import MathJax
    sage: mj = MathJax()
    sage: var('z')
    z
    sage: mj(z^12)
    <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}z^{12}</script></html>
    sage: mj(QQ)
    <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}</script></html>
    sage: mj(ZZ[x])
    <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}[x]</script></html>
    sage: mj(integrate(z^4, z))
    <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}\frac{1}{5} \, z^{5}</script></html>

Basic Use
=========

As indicated in the overview, the simplest way to exploit Sage's
support of LaTeX is to use the ``latex()`` function to create
legitimate LaTeX code to represent mathematical objects.  These
strings can then be incorporated into standalone LaTeX documents.
This works identically in the notebook and at the Sage command
line.

At the other extreme is the ``view()`` command.  At the Sage
command line the command ``view(foo)`` will create the LaTeX
representation of ``foo``, incorporate this into a simple LaTeX
document, and then process that document with your system-wide
TeX installation.  Finally, the appropriate viewer will be called
to display the output from the TeX command.  Which version of TeX
is used, and therefore the nature of the output and associated
viewer, can be customized (see :ref:`sec-custom-processing`).

In the notebook, the ``view(foo)`` command creates the
appropriate combination of HTML and CSS so that MathJax will
render the LaTeX representation properly in the worksheet.  To the
user, it simply creates a nicely formatted version of the output,
distinct from the default ASCII output of Sage.  Not every
mathematical object in Sage has a LaTeX representation amenable to
the limited capabilities of MathJax.  In these cases, the MathJax
interpretation can be bypassed, the system-wide TeX called
instead, and the subsequent output converted to a graphic image
for display in the worksheet.  Affecting and controlling this
process is discussed below in the section
:ref:`sec-custom-generation`.

The internal ``pretty_print()`` command illustrates the
conversion of Sage objects to HTML code that employs MathJax in
the notebook.  ::

    sage: from sage.misc.latex import pretty_print
    sage: pretty_print(x^12)
    <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}x^{12}</script></html>
    sage: pretty_print(integrate(sin(x), x))
    <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}-\cos\left(x\right)</script></html>

The notebook has two other features for employing TeX.
The first is the "Typeset" button just above the first cell of a
worksheet, to the right of the four drop-down boxes.  When
checked, any subsequent evaluations of cells will result in
output interpreted by MathJax, hence of a typeset quality.  Note
that this effect is not retroactive -- previously evaluated cells
need to be re-evaluated.  Essentially, checking the "Typeset"
button is identical to wrapping the output of each cell in the
``view()`` command.

A second feature of the notebook is entering TeX as
part of annotating a worksheet.  When the cursor is placed
between cells of a worksheet so that a blue bar appears, then a
shift-click will open a mini-word-processor, TinyMCE.  This
allows for the entry of text, using a WSIWYG editor to create
HTML and CSS command for styled text.  So it is possible to add
formatted text as commentary within a worksheet.  However, text
between pairs of dollar signs, or pairs of double dollar signs is
interpreted by MathJax as inline or display math (respectively).

.. _sec-custom-generation:

Customizing LaTeX Generation
============================

There are several ways to customize the actual LaTeX code generated by
the ``latex()`` command.  In the notebook and at the Sage command-line
there is a pre-defined object named ``latex`` which has several methods,
which you can list by typing ``latex.``, followed by the tab key
(note the period).

A good example is the ``latex.matrix_delimiters`` method.  It can be
used to change the notation surrounding a matrix -- large parentheses,
brackets, braces, vertical bars.  No notion of style is enforced,
you can mix and match as you please.  Notice how the backslashes
needed in LaTeX require an extra slash so they are escaped
properly within the Python string.  ::

    sage: A = matrix(ZZ, 2, 2, range(4))
    sage: latex(A)
    \left(\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right)
    sage: latex.matrix_delimiters(left='[', right=']')
    sage: latex(A)
    \left[\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right]
    sage: latex.matrix_delimiters(left='\\{', right='\\}')
    sage: latex(A)
    \left\{\begin{array}{rr}
    0 & 1 \\
    2 & 3
    \end{array}\right\}

The ``latex.vector_delimiters`` method works similarly.

The way common rings and fields (integers, rational, reals, etc.)
are typeset can be controlled by the ``latex.blackboard_bold``
method.  These sets are by default typeset in bold, but may
optionally be written in a double-struck fashion as sometimes
done in written work.  This is accomplished by redefining the
``\Bold{}`` macro which is built-in to Sage. ::

    sage: latex(QQ)
    \Bold{Q}
    sage: from sage.misc.latex import MathJax
    sage: mj=MathJax()
    sage: mj(QQ)
    <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Q}</script></html>
    sage: latex.blackboard_bold(True)
    sage: mj(QQ)
    <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbb{#1}}\Bold{Q}</script></html>
    sage: latex.blackboard_bold(False)

It is possible to take advantage of the extensible nature of
TeX by adding in new macros and new packages.  First,
individual macros can be added so that they are used when
MathJax interprets a snippet of TeX in the notebook.  ::

    sage: latex.extra_macros()
    ''
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: latex.extra_macros()
    '\\newcommand{\\foo}{bar}'
    sage: var('x y')
    (x, y)
    sage: latex(x+y)
    x + y
    sage: from sage.misc.latex import MathJax
    sage: mj=MathJax()
    sage: mj(x+y)
    <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}\newcommand{\foo}{bar}x + y</script></html>

Additional macros added this way will also be used in the event
that the system-wide version of TeX is called on
something larger than MathJax can handle.  The command
``latex_extra_preamble`` is used to build the preamble of a
complete LaTeX document, so the following illustrates
how this is accomplished. As usual note the need for the
double-backslashes in the Python strings.  ::


    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: from sage.misc.latex import latex_extra_preamble
    sage: print latex_extra_preamble()
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    sage: latex.add_macro("\\newcommand{\\foo}{bar}")
    sage: print latex_extra_preamble()
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}
    \newcommand{\foo}{bar}

Again, for larger or more complicated LaTeX
expressions, it is possible to add packages (or anything else) to
the preamble of the LaTeX file.  Anything may be
incorporated into the preamble with the ``latex.add_to_preamble``
command, and the specialized command
``latex.add_package_to_preamble_if_available`` will first check
if a certain package is actually available before trying to add
it to the preamble.

Here we add the geometry package to the preamble and use it to
set the size of the region on the page that TeX will
use (effectively setting the margins).  As usual, note the need
for the double-backslashes in the Python strings.  ::


    sage: from sage.misc.latex import latex_extra_preamble
    sage: latex.extra_macros('')
    sage: latex.extra_preamble('')
    sage: latex.add_to_preamble('\\usepackage{geometry}')
    sage: latex.add_to_preamble('\\geometry{letterpaper,total={8in,10in}}')
    sage: latex.extra_preamble()
    '\\usepackage{geometry}\\geometry{letterpaper,total={8in,10in}}'
    sage: print latex_extra_preamble()
    \usepackage{geometry}\geometry{letterpaper,total={8in,10in}}
    \newcommand{\ZZ}{\Bold{Z}}
    ...
    \newcommand{\Bold}[1]{\mathbf{#1}}

A particular package may be added along with a check on its existence,
as follows.  As an example, we just illustrate an attempt to add to
the preamble a package that presumably does not exist. ::

    sage: latex.extra_preamble('')
    sage: latex.extra_preamble()
    ''
    sage: latex.add_to_preamble('\\usepackage{foo-bar-unchecked}')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'
    sage: latex.add_package_to_preamble_if_available('foo-bar-checked')
    sage: latex.extra_preamble()
    '\\usepackage{foo-bar-unchecked}'

.. _sec-custom-processing:

Customizing LaTeX Processing
============================

It is also possible to control which variant of TeX is
used for system-wide invocations, thus also influencing the
nature of the output.  Similarly, it is also possible to control
when the notebook will use MathJax (simple TeX snippets)
or the system-wide TeX installation (more complicated
LaTeX expressions).

The ``latex.engine()`` command can be used to control if the
system-wide executables ``latex``, ``pdflatex`` or ``xelatex``
are employed for more complicated LaTeX expressions.
When ``view()`` is called from the sage command-line and the
engine is set to ``latex``, a dvi file is produced and Sage will
use a dvi viewer (like xdvi) to display the result.  In contrast,
using ``view()`` at the Sage command-line, when the engine is set
to ``pdflatex``, will produce a PDF as the result and Sage will
call your system's utility for displaying PDF files (acrobat,
okular, evince, etc.).

In the notebook, it is necessary to intervene in the decision as
to whether MathJax will interpret a snippet of TeX, or
if the LaTeX is complicated enough that the system-wide
installation of TeX should do the work instead.  The
device is a list of strings, which if any one is discovered in a
piece of LaTeX code signal the notebook to bypass
MathJax and invoke latex (or whichever executable is set by the
``latex.engine()`` command).  This list is managed by the
``latex.add_to_mathjax_avoid_list`` and
``latex.mathjax_avoid_list`` commands. ::

    sage: latex.mathjax_avoid_list([])
    sage: latex.mathjax_avoid_list()
    []
    sage: latex.mathjax_avoid_list(['foo', 'bar'])
    sage: latex.mathjax_avoid_list()
    ['foo', 'bar']
    sage: latex.add_to_mathjax_avoid_list('tikzpicture')
    sage: latex.mathjax_avoid_list()
    ['foo', 'bar', 'tikzpicture']
    sage: latex.mathjax_avoid_list([])
    sage: latex.mathjax_avoid_list()
    []

Suppose a LaTeX expression is produced in the notebook
with ``view()`` or while the "Typeset" button is checked, and
then recognized as requiring the external LaTeX
installation through the "mathjax avoid list."  Then the selected
executable (as specified by ``latex.engine()``) will process the
LaTeX.  However, instead of then spawning an external
viewer (which is the command-line behavior), Sage will attempt to
convert the result into a single, tightly-cropped image, which is
then inserted into the worksheet as the output of the cell.

Just how this conversion proceeds depends on several factors --
mostly which executable you have specified as the engine and
which conversion utilities are available on your system.  Four
useful converters that will cover all eventualities are
``dvips``, ``ps2pdf``, ``dvipng`` and from the ``ImageMagick`` suite,
``convert``.  The goal is to produce a PNG file as the output for
inclusion back into the worksheet.  When a LaTeX
expression can be converted successfully to a dvi by the latex
engine, then dvipng should accomplish the conversion.  If the
LaTeX expression and chosen engine creates a dvi with
specials that dvipng cannot handle, then dvips will create a
PostScript file. Such a PostScript file, or a PDF file created by
an engine such as ``pdflatex``, is then processed into a PNG with
the ``convert`` utility.  The presence of two of these converters
can be tested with the ``have_dvipng()`` and ``have_convert()``
routines.

These conversions are done automatically if you have the necessary
converters installed; if not, then an error message is printed telling
you what's missing and where to download it.

For a concrete example of how complicated LaTeX
expressions can be processed, see the example in the next section
(:ref:`sec-tkz-graph`) for using the LaTeX
``tkz-graph`` package to produce high-quality renderings of
combinatorial graphs.  For other examples, there are some
pre-packaged test cases.  To use these, it is necessary to import
the ``sage.misc.latex.latex_examples`` object, which is an
instance of the ``sage.misc.latex.LatexExamples`` class, as
illustrated below.  This class currently has examples of
commutative diagrams, combinatorial graphs, knot theory and
pstricks, which respectively exercise the following packages:
xy, tkz-graph, xypic, pstricks.  After the import, use
tab-completion on ``latex_examples`` to see the pre-packaged
examples.  Calling each example will give you back some
explanation about what is required to make the example render
properly.  To actually see the examples, it is necessary to use
``view()`` (once the preamble, engine, etc are all set properly).
::

    sage: from sage.misc.latex import latex_examples
    sage: latex_examples.diagram()
    LaTeX example for testing display of a commutative diagram produced
    by xypic.
    <BLANKLINE>
    To use, try to view this object -- it won't work.  Now try
    'latex.add_to_preamble("\\usepackage[matrix,arrow,curve,cmtip]{xy}")',
    and try viewing again -- it should work in the command line but not
    from the notebook.  In the notebook, run
    'latex.add_to_mathjax_avoid_list("xymatrix")' and try again -- you
    should get a picture (a part of the diagram arising from a filtered
    chain complex).

.. _sec-tkz-graph:

An Example: Combinatorial Graphs with tkz-graph
===============================================

High-quality illustrations of combinatorial graphs (henceforth
just "graphs") are possible with the ``tkz-graph`` package.
This package is built on top of the ``tikz`` front-end to the
``pgf`` library.  So all of these components need to be part
of a system-wide TeX installation, and it may be possible
that these components may not be at their most current
versions as packaged in some TeX implementations. So for
best results, it could be necessary or advisable to install
these as part of your personal texmf tree.  Creating,
maintaining and customizing a system-wide or personal TeX
installation is beyond the scope of this document, but it should
be easy to find instructions.  The necessary files are listed in
:ref:`sec-system-wide-tex`.

Thus, to start we need to insure that the relevant packages
are included by adding them to the preamble of the eventual
LaTeX document.  The images of graphs do not form properly
when a dvi file is used as an intermediate format, so it is
best to set the latex engine to the ``pdflatex`` executable.
At this point a command like ``view(graphs.CompleteGraph(4))``
should succeed at the Sage command-line and produce a PDF
with an appropriate image of the complete graph `K_4`.

For a similar experience in the notebook, it is necessary
to disable MathJax processing of the LaTeX code for the graph
by using the "mathjax avoid list."  Graphs are included with a
``tikzpicture`` environment, so this is a good choice for
a string to include in the avoidance list.  Now,
``view(graphs.CompleteGraph(4))`` in a worksheet
should call pdflatex to create a PDF and then the
``convert`` utility will extract a PNG graphic to
insert into the output cell of the worksheet.
The following commands illustrate the steps to get
graphs processed by LaTeX in the notebook. ::

    sage: from sage.graphs.graph_latex import setup_latex_preamble
    sage: setup_latex_preamble()
    sage: latex.extra_preamble() # random - depends on system's TeX installation
    '\\usepackage{tikz}\n\\usepackage{tkz-graph}\n\\usepackage{tkz-berge}\n'
    sage: latex.engine('pdflatex')
    sage: latex.add_to_mathjax_avoid_list('tikzpicture')
    sage: latex.mathjax_avoid_list()
    ['tikzpicture']

At this point, a command like ``view(graphs.CompleteGraph(4))``
should produce a graphic version of the graph pasted into the
notebook, having used ``pdflatex`` to process ``tkz-graph``
commands to realize the graph. Note that there is a variety of
options to affect how a graph is rendered in LaTeX via
``tkz-graph``, which is again outside the scope of this section,
see the section of the Reference manual titled "LaTeX Options for
Graphs" for instructions and details.

.. _sec-system-wide-tex:

A Fully Capable TeX Installation
================================
Many of the more advanced features of the integration of
TeX with Sage requires a system-wide installation of
TeX.  Many versions of Linux have base TeX
packages based on TeX-live, for OSX there is
TeXshop and for Windows there is MikTeX.
The ``convert`` utility is part of the
`ImageMagick <http://www.imagemagick.org/>`_ suite (which
should be a package or an easy download), and the three
programs ``dvipng``, ``ps2pdf``, and ``dvips`` may be
included with your TeX distribution.  The first two may
also be obtained, respectively, from
http://sourceforge.net/projects/dvipng/ and as part of
`Ghostscript <http://www.ghostscript.com/>`_.

Rendering combinatorial graphs requires a recent version of the
PGF library, and the files ``tkz-graph.sty``, ``tkz-arith.sty``
and perhaps ``tkz-berge.sty``, all from the `Altermundus site
<http://altermundus.com/pages/tkz/graph/>`_.

External Programs
=================

There are three programs available to further integrate
TeX and Sage. The first is sagetex.  A concise
description of sagetex is that it is a collection of
TeX macros that allow a LaTeX document to
include instructions to have Sage compute various objects and/or
format objects using the ``latex()`` support built in to Sage.
So as an intermediate step of compiling a LaTeX
document, all of the computational and LaTeX-formatting
features of Sage can be handled automatically.  As an example, a
mathematics examination can maintain a correct correspondence
between questions and answers by using sagetex to have Sage
compute one from the other.  See :ref:`sec-sagetex` for more
information.


tex2sws begins with a LaTeX document, but defines extra
environments for the placement of Sage code.  When processed with
the right tools, the result is a Sage worksheet, with content
properly formatted for MathJax and the Sage code incorporated as
input cells.  So a textbook or article can be authored in
LaTeX, blocks of Sage code included, and the whole
document can be transformed into a Sage worksheet where the
mathematical text is nicely formatted and the blocks of Sage code
are "live."  Currently in development, see `tex2sws @ BitBucket
<http://bitbucket.org/rbeezer/tex2sws/>`_ for more information.

sws2tex reverses the process by beginning with a Sage worksheet
and converting it to legitimate LaTeX for subsequent
processing with all the tools available for LaTeX
documents.  Currently in development, see `sws2tex @ BitBucket
<http://bitbucket.org/whuss/sws2tex/>`_ for more information.
