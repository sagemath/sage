pygambit: Computation in game theory
====================================

Description
-----------

pygambit is the Python interface of the gambit project, a library of
algorithms for computing Nash equilibria of finite strategic and extensive
form games.

In Sage it backs most of the algorithms offered by
:meth:`sage.game_theory.normal_form_game.NormalFormGame.obtain_nash` --
``'LCP'``, ``'lp'``, ``'enummixed'``, ``'enumpure'``, ``'gnm'``, ``'ipa'``,
``'simpdiv'``, ``'logit'``, ``'liap'`` and ``'enumpoly'``, the last one being
the default for games with more than two players -- as well as the ``'gambit'``
MILP solver of the ``'lp'`` algorithm.  It also underlies
:class:`sage.game_theory.extensive_form_game.ExtensiveFormGame`, which wraps a
``pygambit`` game.  Sage's own ``'enumeration'`` and ``'lrs'`` algorithms do
not need it.

See https://www.gambit-project.org/ and
https://gambitproject.readthedocs.io/

License
-------

-  GPL v2 or later


Upstream Contact
----------------

-  https://github.com/gambitproject/gambit

Dependencies
------------

Python.

The ``'enumpoly'`` algorithm can optionally use PHCpack instead of gambit's
built-in polynomial solver; see
https://homepages.math.uic.edu/~jan/download.html

Special Update/Build Instructions
---------------------------------

Make sure corresponding optional doctests still pass:

   sage -t --long --optional=pygambit,sage --all
