.. _sec-sagetex:

*************
Using SageTeX
*************

The SageTeX package allows you to embed the results of Sage computations into a
LaTeX document. To use it, you will need to "install" it first (see
:ref:`sec-sagetex_install`).

An example
----------

Here is a very brief example of using SageTeX. The full documentation
can be found in :file:`SAGE_ROOT/venv/share/doc/sagetex`,
where ``SAGE_ROOT`` is the directory where your Sage installation is
located. That directory contains the documentation and an example file.
See :file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex` for
some possibly useful Python scripts.

To see how SageTeX works, follow the directions for installing SageTeX (in
:ref:`sec-sagetex_install`) and copy the following text into a file named, say,
``st_example.tex``:

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
complaints, which will include:

.. CODE-BLOCK:: text

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
:file:`SAGE_ROOT/venv/share/doc/sagetex`.

.. _sec-sagetex_install:

Make SageTeX known to TeX
-------------------------

Sage is largely self-contained, but some parts do need some intervention
to work properly. SageTeX is one such part.

The SageTeX package allows one to embed computations and plots from Sage
into a LaTeX document. SageTeX is installed in Sage by default, but to
use SageTeX with your LaTeX documents, you need to make your TeX
installation aware of it before it will work.

The key to this is that TeX needs to be able to find ``sagetex.sty``,
which can be found in
:file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex/`, where
``SAGE_ROOT`` is the directory where you built or installed Sage. If
TeX can find ``sagetex.sty``, then SageTeX will work. There are several
ways to accomplish this.

- The first and simplest way is simply to copy ``sagetex.sty`` into the
  same directory as your LaTeX document. Since the current directory is
  always searched when typesetting a document, this will always work.

  There are a couple small problems with this, however: the first is
  that you will end up with many unnecessary copies of ``sagetex.sty``
  scattered around your computer. The second and more serious problem is
  that if you upgrade Sage and get a new version of SageTeX, the Python
  code and LaTeX code for SageTeX may no longer match, causing errors.

- The second way is to use the ``TEXMFLOCAL`` environment variable. If
  you are using the bash shell, you can do

  .. CODE-BLOCK:: shell-session

      $ export TEXMFLOCAL=SAGE_ROOT/venv/share/texmf
      $ mktexlsr       # update kpathsea ls-R databases

  where ``SAGE_ROOT`` is the location of your Sage installation.
  Thereafter, TeX and friends will find the SageTeX style file.
  If you want to make this change persistent, you can add the 1st of the
  above lines to your ``.bashrc`` file.
  If you are using a different shell, you may have to
  modify the above command to make the environment variable known; see
  your shell's documentation for how to do that.

  If you ever move your Sage installation, or install a new version into
  a new directory, you'll need to update the above command to reflect
  the new value of ``SAGE_ROOT``.

- The third (and best) way to make TeX aware of ``sagetex.sty`` is to
  copy that file into a convenient place in your home directory. In most
  TeX distributions, the ``texmf`` directory in your home directory is
  automatically searched for packages. To find out exactly what this
  directory is, do the following on the command line:

  .. CODE-BLOCK:: shell-session

      $ kpsewhich -var-value=TEXMFHOME

  which will print out a directory, such as ``/home/drake/texmf`` or
  ``/Users/drake/Library/texmf``. Copy the ``tex/`` directory from
  :file:`SAGE_ROOT/venv/share/texmf/` into your home ``texmf`` directory
  with a command like

  .. CODE-BLOCK:: shell-session

      $ cp -R SAGE_ROOT/venv/share/texmf/tex TEXMFHOME

  where ``SAGE_ROOT`` is, as usual, replaced with the location of your
  Sage installation and ``TEXMFHOME`` is the result of the
  ``kpsewhich`` command above.

  If you upgrade Sage and discover that SageTeX no longer works, you can
  simply repeat these steps and the Sage and TeX parts of SageTeX will
  again be synchronized.

.. _sagetex_installation_multiuser:

- For installation on a multiuser system, you just modify the above
  instructions appropriately to copy ``sagetex.sty`` into a systemwide
  TeX directory. Instead of the directory ``TEXMFHOME``, probably the
  best choice is to use the result of

  .. CODE-BLOCK:: shell-session

      $ kpsewhich -var-value=TEXMFLOCAL

  which will likely produce something like ``/usr/local/share/texmf``.
  Copy the ``tex`` directory as above into the ``TEXMFLOCAL``
  directory. Now you need to update TeX's database of packages, which
  you can do simply by running

  .. CODE-BLOCK:: shell-session

      $ texhash TEXMFLOCAL

  as root, replacing ``TEXMFLOCAL`` appropriately. Now all users of your
  system will have access to the LaTeX package, and if they can also run
  Sage, they will be able to use SageTeX.

.. warning::

  it's very important that the file ``sagetex.sty`` that LaTeX uses when
  typesetting your document match the version of SageTeX that Sage is
  using. If you upgrade your Sage installation, you really should delete
  all the old versions of ``sagetex.sty`` floating around.

  Because of this problem, we recommend copying the SageTeX files into
  your home directory's texmf directory (the third method above). Then
  there is only one thing you need to do (copy a directory) when you
  upgrade Sage to insure that SageTeX will work properly.

SageTeX documentation
---------------------

While not strictly part of installation, it bears mentioning here that
the documentation for SageTeX is maintained in
:file:`SAGE_ROOT/venv/share/doc/sagetex/sagetex.pdf`. There is also an
example file in the same directory -- see ``example.tex`` and
``example.pdf``, the pre-built result of typesetting that file with
LaTeX and Sage. You can also get those files from the `SageTeX page <https://github.com/sagemath/sagetex>`_.

SageTeX and TeXLive
-------------------

One potentially confusing issue is that the popular TeX distribution
`TeXLive <http://www.tug.org/texlive/>`_ includes SageTeX. This may
seem nice, but with SageTeX, it's important that the Sage bits and LaTeX
bits be synchronized -- which is a problem in this case, since
TeXLive, as shipped by your OS distro, or package manager, might be out of sync with
TeXLive distribution, and the latter might also be out of sync with
the current SageTeX.

Because of this, it is *strongly recommended* that you always install
the LaTeX part of SageTeX from Sage, as described above. The
instructions above will insure that both halves of SageTeX are
compatible and will work properly.
