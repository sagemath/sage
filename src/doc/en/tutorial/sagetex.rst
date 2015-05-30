.. _sec-sagetex:

*************
Using SageTeX
*************

The SageTeX package allows you to embed the results of Sage computations
into a LaTeX document. It comes standard with Sage. To use it, you will
need to "install" it into your local TeX system; here "install" means
copying a single file. See :ref:`installation` in this tutorial and the
"Make SageTeX known to TeX" section of the `Sage installation guide
<http://sagemath.org/doc/installation/index.html>`_ (`this link
<../installation/index.html>`_ should take you to a local copy of the
installation guide) for more information on doing that.

Here is a very brief example of using SageTeX. The full documentation
can be found in ``SAGE_ROOT/local/share/doc/sagetex``,
where ``SAGE_ROOT`` is the directory where your Sage installation is
located. That directory contains the documentation and an example file.
See ``SAGE_ROOT/local/share/texmf/tex/generic/sagetex`` for
some possibly useful Python scripts.

To see how SageTeX works, follow the directions for installing SageTeX
(in :ref:`installation`) and copy the following text into a file named,
say, ``st_example.tex``:

.. warning::

  The text below will have several errors about unknown control
  sequences if you are viewing this in the "live" help. Use the static
  version to see the correct text.

.. code-block:: latex

    \documentclass{article}
    \usepackage{sagetex}

    \begin{document}

    Using Sage\TeX, one can use Sage to compute things and put them into
    your \LaTeX{} document. For example, there are
    $\sage{number_of_partitions(1269)}$ integer partitions of $1269$.
    You don't need to compute the number yourself, or even cut and paste
    it from somewhere.

    Here's some Sage code:

    \begin{sageblock}
        f(x) = exp(x) * sin(2*x)
    \end{sageblock}

    The second derivative of $f$ is

    \[
      \frac{\mathrm{d}^{2}}{\mathrm{d}x^{2}} \sage{f(x)} =
      \sage{diff(f, x, 2)(x)}.
    \]

    Here's a plot of $f$ from $-1$ to $1$:

    \sageplot{plot(f, -1, 1)}

    \end{document}

Run LaTeX on ``st_example.tex`` as usual. Note that LaTeX will have some
complaints, which will include::

    Package sagetex Warning: Graphics file
    sage-plots-for-st_example.tex/plot-0.eps on page 1 does not exist. Plot
    command is on input line 25.

    Package sagetex Warning: There were undefined Sage formulas and/or
    plots. Run Sage on st_example.sagetex.sage, and then run LaTeX on
    st_example.tex again.

Notice that, in addition to the usual collection of files produced by
LaTeX, there is a file called ``st_example.sagetex.sage``. That is a Sage script
produced when you run LaTeX on ``st_example.tex``. The warning message
told you to run Sage on ``st_example.sagetex.sage``, so take its advice and do
that. It will tell you to run LaTeX on ``st_example.tex`` again, but
before you do that, notice that a new file has been created:
``st_example.sagetex.sout``. That file contains the results of Sage's
computations, in a format that LaTeX can use to insert into your text. A
new directory containing an EPS file of your plot has also been created.
Run LaTeX again and you'll see that everything that Sage computed and
plotted is now included in your document.

The different macros used above should be pretty easy to understand. A
``sageblock`` environment typesets your code verbatim and also executes
the code when you run Sage. When you do ``\sage{foo}``, the result put
into your document is whatever you get from running ``latex(foo)``
inside Sage. Plot commands are a bit more complicated, but in their
simplest form, ``\sageplot{foo}`` inserts the image you get from doing
``foo.save('filename.eps')``.

In general, the mantra is:

    - run LaTeX on your .tex file;
    - run Sage on the generated .sage file;
    - run LaTeX again.

You can omit running Sage if you haven't changed around any Sage
commands in your document.

There's a lot more to SageTeX, and since both Sage and LaTeX are
complex, powerful tools, it's a good idea to read the documentation for
SageTeX, which is in
``SAGE_ROOT/local/share/doc/sagetex``.
