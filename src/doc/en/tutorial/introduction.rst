************
Introduction
************

This tutorial should take at most 3-4 hours to fully
work through. You can read it in HTML or PDF versions, or from the
Sage notebook click ``Help``, then click ``Tutorial`` to interactively
work through the tutorial from within Sage.

Though much of Sage is implemented using Python, no Python
background is needed to read this tutorial. You will want to learn
Python (a very fun language!) at some point, and there are many
excellent free resources for doing so: the Python Beginner's Guide [PyB]_
lists many options.  If you just want to quickly try out Sage, this
tutorial is the place to start. For example:

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

.. _installation:

Installation
============

If you do not have Sage installed on a computer and just
want to try some commands, use it online at http://sagecell.sagemath.org.

See the Sage Installation Guide in the documentation section of the
main Sage webpage [SA]_ for instructions on installing Sage on your
computer. Here we merely make a few comments.


#. The Sage download file comes with "batteries included". In other
   words, although Sage uses Python, IPython, PARI, GAP, Singular,
   Maxima, NTL, GMP, and so on, you do not need to install them
   separately as they are included with the Sage distribution.
   However, to use certain Sage features, e.g., Macaulay or KASH, you
   must have the relevant programs installed on your computer already.

#. The pre-compiled binary version of Sage (found on the Sage web
   site) may be easier and quicker to install than the source code
   version. Just unpack the file and run ``sage``.

#. If you'd like to use the SageTeX package (which allows you to embed
   the results of Sage computations into a LaTeX file), you will need to
   make SageTeX known to your TeX distribution. To do this, see the
   section "Make SageTeX known to TeX" in the `Sage installation guide
   <http://doc.sagemath.org/html/en/>`_ (`this link
   <../installation/index.html>`_ should take you to a local copy of the
   installation guide). It's quite easy; you just need to set an
   environment variable or copy a single file to a directory that TeX
   will search.

   The documentation for using SageTeX is located in
   ``$SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/``, where
   "``$SAGE_ROOT``" refers to the directory where you installed Sage --
   for example, ``/opt/sage-9.6``.

Ways to Use Sage
================

You can use Sage in several ways.


-  **Notebook graphical interface:** run ``sage -n jupyter``; see
   `the Jupyter documentation on-line <https://jupyter-notebook.readthedocs.io/en/latest/notebook.html>`_,

-  **Interactive command line:** see :ref:`chapter-interactive_shell`,

-  **Programs:** By writing interpreted and compiled programs in
   Sage (see :ref:`section-loadattach` and :ref:`section-compile`), and

-  **Scripts:** by writing stand-alone Python scripts that use the Sage
   library (see :ref:`section-standalone`).


Longterm Goals for Sage
=======================

-  **Useful**: Sage's intended audience is mathematics students
   (from high school to graduate school), teachers, and research
   mathematicians. The aim is to provide software that can be used to
   explore and experiment with mathematical constructions in algebra,
   geometry, number theory, calculus, numerical computation, etc. Sage
   helps make it easier to interactively experiment with mathematical
   objects.

-  **Efficient:** Be fast. Sage uses highly-optimized mature software
   like GMP, PARI, GAP, and NTL, and so is very fast at certain
   operations.

-  **Free and open source:** The source code must be freely
   available and readable, so users can understand what the system is
   really doing and more easily extend it. Just as mathematicians gain
   a deeper understanding of a theorem by carefully reading or at
   least skimming the proof, people who do computations should be able
   to understand how the calculations work by reading documented
   source code. If you use Sage to do computations in a paper you publish,
   you can rest assured that your readers will always have free access
   to Sage and all its source code, and you are even allowed to archive and
   re-distribute the version of Sage you used.

-  **Easy to compile:** Sage should be easy to compile from source
   for Linux, OS X and Windows users. This provides more flexibility
   for users to modify the system.

-  **Cooperation:** Provide robust interfaces to most other
   computer algebra systems, including PARI, GAP, Singular, Maxima,
   KASH, Magma, Maple, and Mathematica. Sage is meant to unify and extend
   existing math software.

-  **Well documented:** Tutorial, programming guide, reference
   manual, and how-to, with numerous examples and discussion of
   background mathematics.

-  **Extensible:** Be able to define new data types or derive from
   built-in types, and use code written in a range of languages.

-  **User friendly**: It should be easy to understand what
   functionality is provided for a given object and to view
   documentation and source code. Also attain a high level of user
   support.

