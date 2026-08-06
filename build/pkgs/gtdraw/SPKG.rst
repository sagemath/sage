gtdraw: Draw extensive form game trees with TikZ
================================================

Description
-----------

gtdraw is the game tree drawing tool of the gambit project. It turns an
extensive form game -- either a ``pygambit`` game or a file in gambit's
``.efg`` or gtdraw's own ``.ef`` format -- into publication-quality TikZ code,
and can further compile it to a LaTeX document, a PDF, a PNG or an SVG.

In Sage it backs the ``'gtdraw'`` backend of
:meth:`sage.game_theory.extensive_form_game.ExtensiveFormGame.plot`, which
wraps gtdraw's TikZ code in a
:class:`~sage.misc.latex_standalone.TikzPicture`; Sage rather than gtdraw then
compiles and exports it.

See https://www.gambit-project.org/gtdraw/ and, for the prerequisites listed
below, https://www.gambit-project.org/gtdraw/installation/

License
-------

-  GPL v3


Upstream Contact
----------------

-  https://github.com/gambitproject/gtdraw

Dependencies
------------

pygambit.

Rendering a tree -- which is what happens whenever such a picture is displayed,
in a notebook as well as at the Sage command line -- or exporting it to PDF,
PNG or SVG, additionally requires a LaTeX installation providing ``pdflatex``
and TikZ (for example MacTeX, TeX Live or MiKTeX); on top of that PNG export
needs ImageMagick and SVG export needs ``pdftocairo`` (poppler-utils).
Generating the TikZ source itself needs none of these.

Special Update/Build Instructions
---------------------------------

Make sure corresponding optional doctests still pass:

   sage -t --long --optional=gtdraw,pygambit,sage --all
