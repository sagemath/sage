.. highlight:: ipython

.. _chapter-code-basics:

===================
General Conventions
===================

There are many ways to contribute to Sage, including sharing scripts
and Jupyter notebooks that implement new functionality using Sage,
improving to the Sage library, or to working on the many underlying
libraries distributed with Sage, see :ref:`spkg`.
This guide focuses on editing the Sage library itself.

Sage is not just about gathering together functionality. It is about
providing a clear, systematic and consistent way to access a large
number of algorithms, in a coherent framework that makes sense
mathematically. In the design of Sage, the semantics of objects, the
definitions, etc., are informed by how the corresponding objects are
used in everyday mathematics.

To meet the goal of making Sage easy to read, maintain, and improve,
all Python/Cython code that is included with Sage should adhere to the
style conventions discussed in this chapter.


.. _section-coding-python:

Python code style
=================

Follow the standard Python formatting rules when writing code for
Sage, as explained at the following URLs:

* :pep:`0008`
* :pep:`0257`

In particular,

- Use 4 spaces for indentation levels. Do not use tabs as they can
  result in indentation confusion. Most editors have a feature that
  will insert 4 spaces when the :kbd:`Tab` key is hit. Also, many editors
  will automatically search/replace leading tabs with 4 spaces.

- Whitespace before and after assignment and binary operator of the
  lowest priority in the expression::

      i = i + 1
      c = (a+b) * (a-b)

- No whitespace before or after the ``=`` sign if it is used for
  keyword arguments::

      def complex(real, imag=0.0):
          return magic(r=real, i=imag)

- No whitespace immediately inside parenthesis, brackets, and braces::

       spam(ham[1], {eggs: 2})
       [i^2 for i in range(3)]

- Use all lowercase function names with words separated by
  underscores. For example, you are encouraged to write Python
  functions using the naming convention::

      def set_some_value():
          return 1

  Note, however, that some functions do have uppercase letters where
  it makes sense. For instance, the function for lattice reduction by
  the LLL algorithm is called ``Matrix_integer_dense.LLL``.

- Use CamelCase for class names::

      class SomeValue():
          def __init__(self, x):
          self._x  = 1

  and factory functions that mimic object constructors, for example
  ``PolynomialRing`` or::

       def SomeIdentityValue(x):
           return SomeValue(1)


.. _chapter-directory-structure:

Files and directory structure
=============================

Roughly, the Sage directory tree is laid out like this. Note that we
use ``SAGE_ROOT`` in the following as a shortcut for the name of the
directory containing the Sage sources:

.. CODE-BLOCK:: text

    SAGE_ROOT/
        sage          # the Sage launcher
        Makefile      # top level Makefile
        build/        # Sage's build system
            pkgs/     # install, patch, and metadata from spkgs
        src/
            setup.py
            ...
            sage/            # Sage library
                ext_data/    # extra Sage resources (legacy)
            bin/             # the scripts in local/bin that are tracked
        upstream/            # tarballs of upstream sources
        local/               # installed binaries

Python Sage library code goes into ``src/sage/`` and uses the following
conventions. Directory names may be plural (e.g. ``rings``) and file
names are almost always singular (e.g. ``polynomial_ring.py``). Note
that the file ``polynomial_ring.py`` might still contain definitions
of several different types of polynomial rings.

.. NOTE::

   You are encouraged to include miscellaneous notes, emails, design
   discussions, etc., in your package.  Make these plain text files
   (with extension ``.txt``) in a subdirectory called ``notes``.

If you want to create a new directory (`package
<https://docs.python.org/3/tutorial/modules.html#packages>`_) in the
Sage library :sage_root:`src/sage` (say, ``measure_theory``), that
directory will usually contain an empty file ``__init__.py``, which
marks the directory as an ordinary package (see
:ref:`section_namespace_packages`), and also a file ``all.py``,
listing imports from this package that are user-facing and important
enough to be in the global namespace of Sage at startup.  The file
``all.py`` might look like this::

    from .borel_measure import BorelMeasure
    from .banach_tarski import BanachTarskiParadox

but it is generally better to use the :mod:`~sage.misc.lazy_import`
framework::

    from sage.misc.lazy_import import lazy_import
    lazy_import('sage.measure_theory.borel_measure', 'BorelMeasure')
    lazy_import('sage.measure_theory.banach_tarski', 'BanachTarskiParadox')

Then in the file :sage_root:`src/sage/all.py`, add a line ::

    from sage.measure_theory.all import *

Adding new top-level packages below :mod:`sage` should be done
sparingly.  It is often better to create subpackages of existing
packages.

Non-Python Sage source code and small supporting files can be
included in one of the following places:

- In the directory of the Python code that uses that file.  When the
  Sage library is installed, the file will be installed in the same
  location as the Python code. This is referred to as "package data".

  The preferred way to access the data from Python is using the
  `importlib.resources API
  <https://importlib-resources.readthedocs.io/en/latest/using.html>`_,
  in particular the function :func:`importlib.resources.files`.
  Using it, you can:

  - open a resource for text reading: ``fd = files(package).joinpath(resource).open('rt')``
  - open a resource for binary reading: ``fd = files(package).joinpath(resource).open('rb')``
  - read a resource as text: ``text = files(package).joinpath(resource).read_text()``
  - read a resource as bytes: ``bytes = files(package).joinpath(resource).read_bytes()``
  - open an xz-compressed resource for text reading: ``fd = lzma.open(files(package).joinpath(resource).open('rb'), 'rt')``
  - open an xz-compressed resource for binary reading: ``fd = lzma.open(files(package).joinpath(resource).open('rb'), 'rb')``

  If the file needs to be used outside of Python, then the
  preferred way is using the context manager
  :func:`importlib.resources.as_file`. It should be imported in the
  same way as shown above.

- Older code in the Sage library accesses
  the package data in more direct ways. For example,
  :sage_root:`src/sage/interfaces/maxima.py` uses the file
  :sage_root:`src/sage/interfaces/maxima.lisp` at runtime, so it
  refers to it as::

    os.path.join(os.path.dirname(__file__), 'sage-maxima.lisp')

- In an appropriate subdirectory of :sage_root:`src/sage/ext_data/`.
  (At runtime, it is then available in the directory indicated by
  ``SAGE_EXTCODE``).  For example, if ``file`` is placed in
  :sage_root:`src/sage/ext_data/directory/` it can be accessed with ::

    from sage.env import SAGE_EXTCODE
    file = os.path.join(SAGE_EXTCODE, 'directory', 'file')

  This practice is deprecated, see :issue:`33037`.

In all cases, the files must be listed (explicitly or via wildcards) in
the section ``options.package_data`` of the file
:sage_root:`pkgs/sagemath-standard/setup.cfg.m4` (or the corresponding
file of another distribution).

Large data files should not be added to the Sage source tree. Instead, it
is proposed to do the following:

- create a separate git repository and upload them there [2]_,

- add metadata to the repository that make it a pip-installable
  package (distribution package), as explained for example in the
  `Python Packaging User Guide
  <https://packaging.python.org/en/latest/tutorials/packaging-projects/>`_,

- `upload it to PyPI
  <https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives>`_,

- create metadata in ``SAGE_ROOT/build/pkgs`` that make your new
  pip-installable package known to Sage; see :ref:`chapter-packaging`.

For guiding examples of external repositories that host large data
files, see https://github.com/sagemath/conway-polynomials, and
https://github.com/gmou3/matroid-database.

.. [2]

  It is also suggested that the files are compressed, e.g., through
  the command ``xz -e``. They can then be read via a command such as
  ``lzma.open(file, 'rt')``.


Learn by copy/paste
===================

For all of the conventions discussed here, you can find many examples
in the Sage library.  Browsing through the code is helpful, but so is
searching: the functions ``search_src``, ``search_def``, and
``search_doc`` are worth knowing about.  Briefly, from the "sage:"
prompt, ``search_src(string)`` searches Sage library code for the
string ``string``. The command ``search_def(string)`` does a similar
search, but restricted to function definitions, while
``search_doc(string)`` searches the Sage documentation.  See their
docstrings for more information and more options.


Headings of Sage library code files
===================================

The top of each Sage code file should follow this format::

    r"""
    <Short one-line summary that ends with no period>

    <Paragraph description>

    EXAMPLES::

    <Lots and lots of examples>

    AUTHORS:

    - Your Name (2024-01-13): initial version
    - Alice Liddell (2024-05-31): added a method; cleaned docstrings
    - Full name (YYYY-MM-DD): short description

    """

    # ****************************************************************************
    #       Copyright (C) 2024 Your Name <your email>
    #
    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 2 of the License, or
    # (at your option) any later version.
    #                  https://www.gnu.org/licenses/
    # ****************************************************************************

As an example, see :sage_root:`src/sage/rings/integer.pyx`, which contains the
implementation for `\ZZ`. The names of the people who made major contributions
to the file appear in the ``AUTHORS`` section. You can add your name to the
list if you belong to the people, but refrain from being verbose in the
description. The ``AUTHORS`` section shows very rough overview of the history,
especially if a lot of people have been working on that source file. The
authoritative log for who wrote what is always the git repository (see the
output of ``git blame``).

All code included with Sage must be licensed under the GPLv2+ or a
compatible, that is, less restrictive license (e.g. the BSD license).


.. _section-docstrings:

Documentation strings
=====================

.. _section-docstring-function:

The docstring of a function: content
-------------------------------------

**Every** function must have a docstring that includes the following
information. You can use the existing functions of Sage as templates.

-  A **one-sentence description** of the function.

   It must be followed by a blank line and end in a period. It describes the
   function or method's effect as a command ("Do this", "Return that"), not as
   a description like "Returns the pathname ...".

   For methods of a class, it is recommended to refer to the ``self`` argument
   in a descriptive way, unless this leads to a confusion. For example, if
   ``self`` is an integer, then ``this integer`` or ``the integer`` is more
   descriptive, and it is preferable to write

   .. CODE-BLOCK:: rest

       Return whether this integer is prime.

-  A **longer description**.

   This is optional if the one-sentence description does not need
   more explanations.

   Start with assumptions of the object, if there are any. For example,

   .. CODE-BLOCK:: rest

       The poset is expected to be ranked.

   if the function raises an exception when called on a non-ranked poset.

   Define your terms

   .. CODE-BLOCK:: rest

       The lexicographic product of `G` and `H` is the graph with vertex set ...

   and mention possible aliases

   .. CODE-BLOCK:: rest

       The tensor product is also known as the categorical product and ...

-  An **INPUT** and an **OUTPUT** block describing the input/output of
   the function.

   The INPUT block describes all arguments that the function accepts.

   1. The type names should be descriptive, but do not have to represent the
      exact Sage/Python types. For example, use "integer" for anything that
      behaves like an integer, rather than "int" or "Integer".

   2. Mention the default values of the input arguments when applicable.

   .. CODE-BLOCK:: rest

       INPUT:

       - ``n`` -- integer

       - ``p`` -- prime integer (default: `2`); coprime with `n`

       - ``var`` -- string (default: ``'lambda'``)

       - ``check`` -- boolean (default: ``True``); specifies whether to check for primality

       - ``algorithm`` -- (default: ``None``) the name of the algorithm to use

   The OUTPUT block describes the expected output. This is required if the
   one-sentence description of the function needs more explanation.

   .. CODE-BLOCK:: rest

       OUTPUT: the plaintext decrypted from the ciphertext ``C``

   It is often the case that the output consists of several items.

   .. CODE-BLOCK:: rest

       OUTPUT: a tuple of

       - the reduced echelon form `H` of the matrix `A`

       - the transformation matrix `U` such that `UA = H`

   You are recommended to be verbose enough for complicated outputs.

   .. CODE-BLOCK:: rest

       OUTPUT:

       The decomposition of the free module on which this matrix `A` acts from
       the right (i.e., the action is `x` goes to `xA`), along with whether
       this matrix acts irreducibly on each factor. The factors are guaranteed
       to be sorted in the same way as the corresponding factors of the
       characteristic polynomial.

-  An **EXAMPLES** block for examples. This is not optional.

   These examples are used for documentation, but they are also
   tested before each release just like TESTS block.

   They should have good coverage of the functionality in question.

-  A **SEEALSO** block (highly recommended) with links to related parts of
   Sage. This helps users find the features that interest them and discover
   the new ones.

   .. CODE-BLOCK:: rest

       .. SEEALSO::

           :ref:`chapter-sage_manuals_links`,
           :meth:`sage.somewhere.other_useful_method`,
           :mod:`sage.some.related.module`.

   See :ref:`chapter-sage_manuals_links` for details on how to setup
   links in Sage.

-  An **ALGORITHM** block (optional).

   It indicates what algorithm and/or what software is used, e.g.
   ``ALGORITHM: Uses Pari``. Here's a longer example with a
   bibliographical reference:

   .. CODE-BLOCK:: rest

       ALGORITHM:

       The following algorithm is adapted from page 89 of [Nat2000]_.

       Let `p` be an odd (positive) prime and let `g` be a generator
       modulo `p`. Then `g^k` is a generator modulo `p` if and only if
       `\gcd(k, p-1) = 1`. Since `p` is an odd prime and positive, then
       `p - 1` is even so that any even integer between 1 and `p - 1`,
       inclusive, is not relatively prime to `p - 1`. We have now
       narrowed our search to all odd integers `k` between 1 and `p - 1`,
       inclusive.

       So now start with a generator `g` modulo an odd (positive) prime
       `p`. For any odd integer `k` between 1 and `p - 1`, inclusive,
       `g^k` is a generator modulo `p` if and only if `\gcd(k, p-1) = 1`.

   The bibliographical reference should go in Sage's master
   bibliography file,
   :sage_root:`src/doc/en/reference/references/index.rst`:

   .. CODE-BLOCK:: rest

       .. [Nat2000] \M. B. Nathanson. Elementary Methods in Number Theory.
          Springer, 2000.

-  A **NOTE** block for tips/tricks (optional).

   .. CODE-BLOCK:: rest

       .. NOTE::

           You should note that this sentence is indented at least 4
           spaces. Never use the tab character.

- A **WARNING** block for critical information about your code (optional).

  For example known situations for which the code breaks, or anything
  that the user should be aware of.

  .. CODE-BLOCK:: rest

      .. WARNING::

          Whenever you edit the Sage documentation, make sure that
          the edited version still builds. That is, you need to ensure
          that you can still build the HTML and PDF versions of the
          updated documentation. If the edited documentation fails to
          build, it is very likely that you would be requested to
          change your patch.

- A **TODO** block for future improvements (optional).

  It can contain disabled doctests to demonstrate the desired
  feature. Here's an example of a TODO block:

  .. CODE-BLOCK:: rest

      .. TODO::

          Add to ``have_fresh_beers`` an interface with the faster
          algorithm "Buy a Better Fridge" (BaBF)::

              sage: have_fresh_beers('Bière de l\'Yvette', algorithm="BaBF") # not implemented
              Enjoy !

- A **PLOT** block to illustrate with pictures the output of a function.

  Generate with Sage code an object ``g`` with a ``.plot`` method, then call
  ``sphinx_plot(g)``:

  .. CODE-BLOCK:: rest

      .. PLOT::

          g = graphs.PetersenGraph()
          sphinx_plot(g)

- A **REFERENCES** block to list related books or papers (optional).

  Almost all bibliographic information should be put in the master bibliography
  file, see below. Citations will then link to the master bibliography where
  the reader can find the bibliographic details (see below for citation
  syntax).  REFERENCE blocks in individual docstrings are therefore usually not
  necessary.

  Nevertheless, a REFERENCE block can be useful if there are relevant sources
  which are not explicitly mentioned in the docstring or if the docstring is
  particularly long. In that case, add the bibliographic information to the
  master bibliography file, if not already present, and add a reference block
  to your docstring as follows:

  .. CODE-BLOCK:: rest

      REFERENCES:

      For more information, see [Str1969]_, or one of the following references:

      - [Sto2000]_

      - [Voe2003]_

  Note the trailing underscores which makes the citations into hyperlinks. See
  below for more about the master bibliography file. For more about citations,
  see the `Sphinx/reST markup for citations
  <https://www.sphinx-doc.org/rest.html#citations>`_. For links to GitHub issues and PRs
  or wikipedia, see :ref:`chapter-sage_manuals_links`.

- A **TESTS** block (highly recommended).

  Formatted just like EXAMPLES, containing tests that are not relevant
  to users.  In particular, these blocks are not shown when users ask
  for help via ``foo?``: they are stripped by the function
  :func:`sage.misc.sagedoc.skip_TESTS_block`.

  Special and corner cases, like number zero, one-element group etc.
  should usually go to this block. This is also right place for most
  tests of input validation; for example if the function accepts
  ``direction='up'`` and ``direction='down'``, you can use this block to check
  that ``direction='junk'`` raises an exception.

  For the purposes of removal, A "TESTS" block is a block starting
  with "TESTS:" (or the same with two colons), on a line on
  its own, and ending either with a line indented less than "TESTS",
  or with a line with the same level of indentation -- not more --
  matching one of the following:

  - a Sphinx directive of the form ".. foo:", optionally followed by
    other text.

  - text of the form "UPPERCASE:", optionally followed by other
    text.

  - lines which look like a reST header: one line containing
    anything, followed by a line consisting only of whitespace,
    followed by a string of hyphens, equal signs, or other
    characters which are valid markers for reST
    headers: ``- = ` : ' " ~ _ ^ * + # < >``.
    However, lines only containing double colons `::` do not
    end "TESTS" blocks.

  Sometimes (but rarely) one has private or protected methods that don't need a
  proper ``EXAMPLES`` doctest. In these cases, one can either write traditional
  doctest using the ``TESTS`` block or use pytest to test the method.
  In the latter case, one has to add ``TESTS: pytest`` to the docstring, so that
  the method is explicitly marked as tested.

Note about Sphinx directives vs. other blocks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The main Sphinx directives that are used in Sage are:

``.. MATH::``, ``.. NOTE::``, ``.. PLOT::``, ``.. RUBRIC::``,
``.. SEEALSO::``, ``.. TODO::``, ``.. TOPIC::`` and ``.. WARNING::``.

They must be written exactly as above, so for example
``WARNING::`` or ``.. WARNING ::`` will not work.

Some other directives are also available, but less frequently used, namely:

``.. MODULEAUTHOR::``, ``.. automethod::``, ``.. autofunction::``,
``.. image::``, ``.. figure::``.

Other blocks shall not be used as directives; for example
``.. ALGORITHM::`` will not be shown at all.

Sage documentation style
^^^^^^^^^^^^^^^^^^^^^^^^

All Sage documentation is written in reStructuredText (reST) and is
processed by Sphinx. See https://www.sphinx-doc.org/rest.html for an
introduction. Sage imposes these styles:

- Lines should be shorter than 80 characters. If in doubt, read `PEP8: Maximum
  Line Length <https://www.python.org/dev/peps/pep-0008/#maximum-line-length>`_.

- All reST and Sphinx directives (like ``.. WARNING::``, ``.. NOTE::``,
  ``.. MATH::``, etc.) are written in uppercase.

- Code fragments are quoted with double backticks. This includes function
  arguments and the Python literals like ````True````, ````False```` and
  ````None````. For example:

  .. CODE-BLOCK:: rest

      If ``check`` is ``True``, then ...

Sage's master **BIBLIOGRAPHY** file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All bibliographical references should be stored in the master
bibliography file,
:sage_root:`src/doc/en/reference/references/index.rst`, in the
format

.. CODE-BLOCK:: rest

  .. [Gau1801] \C. F. Gauss, *Disquisitiones Arithmeticae*, 1801.

  .. [RSA1978] \R. Rivest, A. Shamir, L. Adleman,
               "A Method for Obtaining Digital Signatures and
               Public-Key Cryptosystems".
               Communications of the ACM **21** (February 1978),
               120–126. :doi:`10.1145/359340.359342`.

The part in brackets is the citation key: given these examples, you
could then use ``[Gau1801]_`` in a docstring to provide a link to the
first reference. Note the trailing underscore which makes the citation a
hyperlink.

When possible, the key should have this form: for a single author, use the
first three letters of the family name followed by the year; for multiple
authors, use the first letter of each of the family names followed by the
year. Note that the year should be four digits, not just the last two -- Sage
already has references from both 1910 and 2010, for example.

When abbreviating the first name of an author in a bibliography
listing, be sure to put a backslash in front of it. This ensures
that the letter (``C.`` in the example above) will not be
interpreted as a list enumerator.

For more about citations, see the `Sphinx/reST markup for citations
<https://www.sphinx-doc.org/rest.html#citations>`_.

Template
^^^^^^^^

Use the following template when documenting functions. Note the
indentation:

.. skip    # do not doctest

.. CODE-BLOCK:: python

    def point(self, x=1, y=2):
        r"""
        Return the point `(x^5, y)`.

        INPUT:

        - ``x`` -- integer (default: `1`); the description of the
          argument ``x`` goes here. If it contains multiple lines, all
          the lines after the first need to begin at the same indentation
          as the backtick.

        - ``y`` -- integer (default: `2`); the description of the
          argument ``y``

        OUTPUT: tuple; further description of the output

        EXAMPLES:

        This example illustrates ... ::

            sage: A = EuclideanSpace(2)
            sage: A.point(2, 3)
            (2, 3)

        We now ... ::

            sage: B = A.point(5, 6)
            sage: ...

        It is an error to ... ::

            sage: C = A.point('x', 7)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'x' to an integer

        .. NOTE::

            This function uses :func:`pow` to determine the fifth
            power of `x`.

        ...

        .. SEEALSO::

            :func:`line`

        TESTS::

            sage: A.point(42, 0)  # check for corner case y = 0
            ...
        """
        <body of the function>

The master bibliography file would contain

.. CODE-BLOCK:: rest

        .. [BCDT2001] Breuil, Conrad, Diamond, Taylor,
                      "Modularity ...."

You are strongly encouraged to:

- Use LaTeX typesetting (see :ref:`section-latex-typeset`).

- Use raw strings (``r"""..."""``), regardless of whether the docstring
  currently contains any backslashes or not.

- Liberally describe what the examples do.

  .. NOTE::

     There must be a blank line after the example code and before the
     explanatory text for the next example (indentation is not enough).

- Illustrate the exceptions raised by the function with examples (as
  given above: "It is an error to [..]", ...)

- Include many examples.

  They are helpful for the users, and are crucial for the quality and
  adaptability of Sage. Without such examples, small changes to one part
  of Sage that break something else might not go seen until much later
  when someone uses the system, which is unacceptable.

Fine points on styles
---------------------

A Sage developer, in writing code and docstrings, should follow the styles
suggested in this manual, except special cases with good reasons. However, there
are some details where we as a community did not reach to an agreement on
the official style. These are

- one space::

    This is the first sentence. This is the second sentence.

  vs two spaces::

    This is the first sentence.  This is the second sentence.

  between sentences.

- tight list::

    - first item
    - second item
    - third item

  vs spaced list::

    - first item

    - second item

    - third item

There are different opinions on each of these, and in reality, we find
instances in each style in our codebase. Then what should we do? Do we decide
on one style by voting? There are different opinions even on what to do!

We can at least do this to prevent any dispute about these style conflicts:

- Acknowledge different authors may have different preferences on these.

- Respect the style choice of the author who first wrote the code or the docstrings.

Private functions
^^^^^^^^^^^^^^^^^

Functions whose names start with an underscore are considered
private. They do not appear in the reference manual, and their docstring
should not contain any information that is crucial for Sage users. You
can make their docstrings be part of the documentation of another
method. For example::

    class Foo(SageObject):

        def f(self):
            """
            <usual docstring>

            .. automethod:: _f
            """
            return self._f()

        def _f(self):
             """
             This would be hidden without the ``.. automethod::``
             """

Private functions should contain an EXAMPLES (or TESTS) block.

A special case is the constructor ``__init__``: due to its special
status the ``__init__`` docstring is used as the class docstring if
there is not one already. That is, you can do the following:

.. CODE-BLOCK:: ipycon

    sage: class Foo(SageObject):
    ....:     # no class docstring
    ....:     def __init__(self):
    ....:         """Construct a Foo."""
    sage: foo = Foo()
    sage: from sage.misc.sageinspect import sage_getdoc
    sage: sage_getdoc(foo)              # class docstring
    'Construct a Foo.\n'
    sage: sage_getdoc(foo.__init__)     # constructor docstring
    'Construct a Foo.\n'

.. _section-latex-typeset:

LaTeX typesetting
-----------------

In Sage's documentation LaTeX code is allowed and is marked with **backticks**:

    ```x^2 + y^2 = 1``` yields `x^2 + y^2 = 1`.

**Backslashes:** For LaTeX commands containing backslashes, either use double
backslashes or begin the docstring with a ``r"""`` instead of ``"""``::

    def cos(x):
        """
        Return `\\cos(x)`.
        """

    def sin(x):
        r"""
        Return `\sin(x)`.
        """

We strongly suggest to use the latter.

**MATH block:** This is similar to the LaTeX syntax ``\[<math expression>\]``
(or ``$$<math expression>$$``). For instance:

.. CODE-BLOCK:: rest

    .. MATH::

        \sum_{i=1}^{\infty} (a_1 a_2 \cdots a_i)^{1/i}
        \leq
        e \sum_{i=1}^{\infty} a_i

.. MATH::

    \sum_{i=1}^{\infty} (a_1 a_2 \cdots a_i)^{1/i}
    \leq
    e \sum_{i=1}^{\infty} a_i

The **aligned** environment works as it does in LaTeX:

.. CODE-BLOCK:: rest

    .. MATH::

        \begin{aligned}
         f(x) & = x^2 - 1 \\
         g(x) & = x^x - f(x - 2)
        \end{aligned}

.. MATH::

    \begin{aligned}
     f(x) & = x^2 - 1 \\
     g(x) & = x^x - f(x - 2)
    \end{aligned}

When building the PDF documentation, everything is translated to LaTeX
and each MATH block is automatically wrapped in a math environment --
in particular, it is turned into ``\begin{gather} block
\end{gather}``.  So if you want to use a LaTeX environment (like
``align``) which in ordinary LaTeX would not be wrapped like this, you
must add a **:nowrap:** flag to the MATH mode. See also `Sphinx's
documentation for math blocks
<http://sphinx-doc.org/latest/ext/math.html?highlight=nowrap#directive-math>`_. :

.. CODE-BLOCK:: rest

    .. MATH::
       :nowrap:

       \begin{align}
          1+...+n &= n(n+1)/2\\
          &= O(n^2)\\
       \end{align}

.. MATH::
   :nowrap:

   \begin{align}
   1+...+n &= n(n+1)/2\\
   &= O(n^2)\\
   \end{align}

**Readability balance:** in the interactive console, LaTeX formulas contained
in the documentation are represented by their LaTeX code (with backslashes
stripped). In this situation ``\\frac{a}{b}`` is less readable than ``a/b`` or
``a b^{-1}`` (some users may not even know LaTeX code). Make it pleasant for
everybody as much as you can manage.

**Commons rings** `(\Bold{Z},\Bold{N},...)`: The Sage LaTeX style is to typeset
standard rings and fields using the locally-defined macro ``\\Bold`` (e.g.
``\\Bold{Z}`` gives `\Bold{Z}`).

**Shortcuts** are available which preserve readability, e.g. ``\\ZZ`` (`\ZZ`),
``\\RR`` (`\RR`), ``\\CC`` (`\CC`), and ``\\QQ`` (`\QQ`). They appear as
LaTeX-formatted ``\\Bold{Z}`` in the html manual, and as ``Z`` in the
interactive help. Other examples: ``\\GF{q}`` (`\GF{q}`) and ``\\Zmod{p}``
(`\Zmod{p}`).

See the file :sage_root:`src/sage/misc/latex_macros.py` for a full list and
for details about how to add more macros.

.. _section-doctest-writing:

Writing testable examples
-------------------------

The examples from Sage's documentation have a double purpose:

- They provide **illustrations** of the code's usage to the users

- They are **tests** that are checked before each release, helping us avoid
  new bugs.

All new doctests added to Sage should **pass all tests** (see
:ref:`chapter-doctesting`), i.e. running ``sage -t your_file.py`` should not
give any error messages. Below are instructions about how doctests should be
written.

.. highlight:: ipycon

**What doctests should test:**

- **Interesting examples** of what the function can do. This will be the
  most helpful to a lost user. It is also the occasion to check famous
  theorems (just in case)::

    sage: is_prime(6) # 6 is not prime
    False
    sage: 2 * 3 # and here is a proof
    6

- All **meaningful combinations** of input arguments. For example a function
  may accept an ``algorithm="B"`` argument, and doctests should involve both
  ``algorithm="A"`` and ``algorithm="B"``.

- **Corner cases:** the code should be able to handle a 0 input, or an empty
  set, or a null matrix, or a null function, ... All corner cases should be
  checked, as they are the most likely to be broken, now or in the future. This
  probably belongs to the TESTS block (see :ref:`section-docstring-function`).

- **Interruption:** if the function might take a very long time, use
  :func:`~sage.doctest.util.ensure_interruptible_after` to check that the user
  can interrupt it. For example, the following tests ``sleep(3)`` can be
  interrupted after 1 second::

    sage: from sage.doctest.util import ensure_interruptible_after
    sage: with ensure_interruptible_after(1) as data: sleep(3)

- **Systematic tests** of all small-sized inputs, or tests of **random**
  instances if possible.

  .. NOTE::

     Note that **TestSuites** are an automatic way to generate some of these
     tests in specific situations. See
     :sage_root:`src/sage/misc/sage_unittest.py`.

**The syntax:**

- **Environment:** doctests should work if you copy/paste them in Sage's
  interactive console. For example, the function ``AA()`` in the file
  :sage_root:`src/sage/algebras/steenrod/steenrod_algebra.py` includes an
  EXAMPLES block containing the following::

    sage: from sage.algebras.steenrod.steenrod_algebra import AA as A
    sage: A()
    mod 2 Steenrod algebra, milnor basis

  Sage does not know about the function ``AA()`` by default, so it needs to be
  imported before it is tested. Hence the first line in the example.

  All blocks within the same docstring are linked: Variables set
  in a doctest keep their values for the remaining doctests within the
  same docstring. It is good practice to use different variable names for different
  values, as it makes the data flow in the examples easier to understand
  for human readers.  (It also makes the data flow analysis in the
  Sage doctester more precise.)  In particular, when unrelated examples
  appear in the same docstring, do not use the same variable name
  for both examples.

- **Preparsing:** As in Sage's console, `4/3` returns `4/3` and not
  `1.3333333333333333` as in Python. Testing occurs with full Sage
  preparsing of input within the standard Sage shell environment, as
  described in :ref:`section-preparsing`.

- **Writing files:** If a test outputs to a file, the file should be a
  temporary file.  Use :func:`tmp_filename` to get a temporary filename, or
  :func:`tmp_dir` to get a temporary directory. An example from
  :sage_root:`src/sage/plot/graphics.py`)::

      sage: plot(x^2 - 5, (x, 0, 5), ymin=0).save(tmp_filename(ext='.png'))

- **Multiline doctests:** You may write tests that span multiple lines, using
  the line continuation marker ``....:`` ::

      sage: for n in srange(1,10):
      ....:     if n.is_prime():
      ....:         print(n)
      2
      3
      5
      7

- **Wrap long doctest lines:** Note that all doctests in EXAMPLES blocks
  get formatted as part of our HTML and PDF reference manuals. Our HTML manuals
  are formatted using the responsive design provided by the
  :ref:`Furo theme <spkg_furo>`. Even when the browser window is expanded to
  make use of the full width of a wide desktop screen, the style will not
  allow code boxes to grow arbitrarily wide.

  It is best to wrap long lines when possible so that readers do not have to
  scroll horizontally (back and forth) to follow an example.

  - Try to wrap long lines somewhere around columns 80 to 88
    and try to never exceed column 95 in the source file.
    (Columns numbers are from the left margin in the source file;
    these rules work no matter how deep the docstring may be nested
    because also the formatted output will be nested.)

  - If you have to break an expression at a place that is not already
    nested in parentheses, wrap it in parentheses::

      sage: (len(list(Permutations(['a', 'b', 'c', 'd', 'e', 'f', 'g'])))
      ....:    == len(list(Permutations(7))))
      True

  - If the output in your only example is very wide and cannot be reasonably
    reformatted to fit (for example, large symbolic matrices or numbers with many digits),
    consider showing a smaller example first.

  - No need to wrap long ``import`` statements. Typically, the ``import`` statements
    are not the interesting parts of the doctests. Users only need to be able to
    copy-paste them into a Sage session or source file::

      sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict, MPolynomialRing_polydict_domain  # this is fine

  - Wrap and indent long output to maximize readability in the source code
    and in the HTML output. But do not wrap strings::

      sage: from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_quasi
      sage: P.<x, y, z> = ProjectiveSpace(2, ZZ)
      sage: S = P.subscheme([])
      sage: T = P.subscheme([x - y])
      sage: U = AlgebraicScheme_quasi(S, T); U
      Quasi-projective subscheme X - Y of Projective Space of dimension 2
       over Integer Ring,
        where X is defined by: (no polynomials)
          and Y is defined by: x - y
      sage: U._repr_()                                                                                                                                                    # this is fine
      'Quasi-projective subscheme X - Y of Projective Space of dimension 2 over Integer Ring, where X is defined by:\n  (no polynomials)\nand Y is defined by:\n  x - y'

    Also, if there is no whitespace in the doctest output where you could wrap the line,
    do not add such whitespace. Just don't wrap the line::

      sage: B47 = RibbonGraph(4,7, bipartite=True); B47
      Ribbon graph of genus 9 and 1 boundary components
      sage: B47.sigma()                                                                                                                                                           # this is fine
      (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)(15,16,17,18,19,20,21)(22,23,24,25,26,27,28)(29,30,31,32)(33,34,35,36)(37,38,39,40)(41,42,43,44)(45,46,47,48)(49,50,51,52)(53,54,55,56)

  - Doctest tags for modularization purposes such as ``# needs sage.modules``
    (see :ref:`section-further_conventions`) should be aligned at column 88.
    Clean lines from consistent alignment help reduce visual clutter.
    Moreover, at the maximum window width, only the word ``# needs`` will be
    visible in the HTML output without horizontal scrolling, striking a
    thoughtfully chosen balance between presenting
    the information and reducing visual clutter. (How much can be seen may be
    browser-dependent, of course.) In visually dense doctests, you can try to sculpt out visual space to separate
    the test commands from the annotation.

  - Doctest tags such as ``# optional - pynormaliz`` that make the doctest
    conditional on the presence of optional packages, on the other hand,
    should be aligned so that they are visible without having to scroll horizontally.
    The :ref:`doctest fixer <section-fixdoctests-optional-needs>` uses
    tab stops at columns 48, 56, 64, ... for these tags.

- **Split long lines:** You may want to split long lines of code with a
  backslash. Note: this syntax is non-standard and may be removed in the
  future::

      sage: n = 123456789123456789123456789\
      ....:     123456789123456789123456789
      sage: n.is_prime()
      False

- **Doctests flags:** flags are available to change the behaviour of doctests:
  see :ref:`section-further_conventions`.

.. _section-further_conventions:

Special markup to influence doctests
------------------------------------

Overly complicated output in the example code can be shortened
by an ellipsis marker ``...``::

    sage: [ZZ(n).ordinal_str() for n in range(25)]
    ['0th',
     '1st',
     '2nd',
     '3rd',
     '4th',
     '5th',
     ...
     '21st',
     '22nd',
     '23rd',
     '24th']
    sage: ZZ('sage')
    Traceback (most recent call last):
    ...
    TypeError: unable to convert 'sage' to an integer

On the proper usage of the ellipsis marker, see :python:`Python's documentation
<library/doctest.html#doctest.ELLIPSIS>`.

There are a number of magic comments that you can put into the example
code that change how the output is verified by the Sage doctest
framework. Here is a comprehensive list:

- **random:** The line will be executed, but its output will not be checked with
  the output in the documentation string::

      sage: c = CombinatorialObject([1,2,3])
      sage: hash(c)  # random
      1335416675971793195
      sage: hash(c)  # random
      This doctest passes too, as the output is not checked

  Doctests are expected to pass with any state of the pseudorandom number
  generators (PRNGs).
  When possible, avoid the problem, e.g.: rather than checking the value of the
  hash in a doctest, one could illustrate successfully using it as a key in a
  dict.

  One can also avoid the ``random``-tag by checking basic properties::

      sage: QQ.random_element().parent() is QQ
      True
      sage: QQ.random_element() in QQ
      True
      sage: a = QQ.random_element()
      sage: b = QQ._random_nonzero_element()
      sage: c = QQ._random_nonzero_element()
      sage: (a/c) / (b/c) == a/b
      True

  Distribution can be checked with loops::

      sage: found = {i: False for i in range(-2, 3)}
      sage: while not all(found.values()):
      ....:     found[ZZ.random_element(-2, 3)] = True

  This is mathematically correct, as it is
  guaranteed to terminate. However, there is a
  nonzero probability of a timeout.

- **long time:** The line is only tested if the ``--long`` option is given, e.g.
  ``sage -t --long f.py``.

  Use it for doctests that take more than a second to run. No example should
  take more than about 30 seconds::

      sage: E = EllipticCurve([0, 0, 1, -1, 0])
      sage: E.regulator()        # long time (1 second)
      0.0511114082399688

- **tol** or **tolerance:** The numerical values returned by the line are only
  verified to the given tolerance. It is useful when the output is subject to
  numerical noise due to system-dependent (floating point arithmetic, math
  libraries, ...) or non-deterministic algorithms.

  - This may be prefixed by ``abs[olute]`` or ``rel[ative]`` to specify whether
    to measure **absolute** or **relative** error (see the
    :wikipedia:`Approximation_error`).

  - If none of ``abs/rel`` is specified, the error is considered to be
    ``absolute`` when the expected value is **zero**, and is ``relative`` for
    **nonzero** values.

  ::

     sage: n(pi)  # abs tol 1e-9
     3.14159265358979
     sage: n(pi)  # rel tol 2
     6
     sage: n(pi)  # abs tol 1.41593
     2
     sage: K.<zeta8> = CyclotomicField(8)
     sage: N(zeta8)  # absolute tolerance 1e-10
     0.7071067812 + 0.7071067812*I

  **Multiple numerical values:** the representation of complex numbers,
  matrices, or polynomials usually involves several numerical values. If a
  doctest with tolerance contains several numbers, each of them is checked
  individually::

      sage: print("The sum of 1 and 1 equals 5")  # abs tol 1
      The sum of 2 and 2 equals 4
      sage: e^(i*pi/4).n()  # rel tol 1e-1
      0.7 + 0.7*I
      sage: ((x+1.001)^4).expand()  # rel tol 2
      x^4 + 4*x^3 + 6*x^2 + 4*x + 1
      sage: M = matrix.identity(3) + random_matrix(RR,3,3)/10^3
      sage: M^2 # abs tol 1e-2
      [1 0 0]
      [0 1 0]
      [0 0 1]

  The values that the doctesting framework involves in the error computations
  are defined by the regular expression ``float_regex`` in
  :mod:`sage.doctest.parsing`.

- **not implemented** or **not tested:** The line is never tested.

  Use it for very long doctests that are only meant as documentation. It can
  also be used for todo notes of what will eventually be implemented::

      sage: factor(x*y - x*z)    # not implemented

  It is also immediately clear to the user that the indicated example
  does not currently work.

  .. NOTE::

     Skip all doctests of a file/directory

     - **file:** If one of the first 10 lines of a file starts with any of
       ``r""" nodoctest`` (or ``""" nodoctest`` or ``# nodoctest`` or ``%
       nodoctest`` or ``.. nodoctest``, or any of these with different spacing),
       then that file will be skipped.

     - **directory:** If a directory contains a file ``nodoctest.py``, then that
       whole directory will be skipped.

     Neither of this applies to files or directories which are explicitly given
     as command line arguments: those are always tested.

- **optional** or **needs:** A line tagged with ``optional - FEATURE`` or
  ``needs FEATURE`` is tested if the feature is available in Sage. If
  ``FEATURE`` starts with an exclamation point ``!``, then the condition is
  negated, that is, the doctest runs only if the feature is not available.

  If the feature is included in the ``--optional=KEYWORD`` flag passed to
  ``sage -t`` (see :ref:`section-optional-doctest-flag`), then the line is
  tested regardless of the feature availability.

  The main applications are:

  - **optional packages:** When a line requires an optional package to be
    installed (e.g. the ``rubiks`` package)::

      sage: C = RubiksCube("R*L")
      sage: C.solve()                    # optional - rubiks (a hybrid algorithm is used)
      'L R'
      sage: C.solve()                    # optional - !rubiks (GAP is used)
      'L*R'

  - **features:** When a line requires a feature to be present::

      sage: SloaneEncyclopedia[60843]    # optional - sloane_database
      [1, 6, 21, 107, 47176870]

      sage: SloaneEncyclopedia[60843]    # optional - !sloane_database
      Traceback (most recent call last):
      ...
      OSError: The Sloane Encyclopedia database must be installed. Use e.g.
      'SloaneEncyclopedia.install()' to download and install it.

    For lines that require an internet connection::

       sage: oeis(60843)                 # optional - internet
       A060843: Busy Beaver problem: a(n) = maximal number of steps that an
       n-state Turing machine can make on an initially blank tape before
       eventually halting.

  - **known bugs:** For lines that describe known bugs, you can use ``# optional - bug``,
    although ``# known bug`` is preferred.

    .. CODE-BLOCK:: rest

        The following should yield 4.  See :issue:`2`. ::

            sage: 2+2  # optional - bug
            5
            sage: 2+2  # known bug
            5

  - **modularization:** To enable
    :ref:`separate testing of the distribution packages <section-doctesting-venv>`
    of the modularized Sage library, doctests that depend on features provided
    by other distribution packages can be tagged ``# needs FEATURE``.
    For example:

    .. CODE-BLOCK:: rest

        Consider the following calculation::

            sage: a = AA(2).sqrt()  # needs sage.rings.number_field
            sage: b = sqrt(3)       # needs sage.symbolic
            sage: a + AA(b)         # needs sage.rings.number_field sage.symbolic
            3.146264369941973?

  .. NOTE::

      - Any words after ``# optional`` and ``# needs``  are interpreted as a list of
        package (spkg) names or other feature tags, separated by spaces.

      - Any punctuation other than underscores (``_``) and periods (``.``),
        that is, commas, hyphens, semicolons, ..., after the
        first word ends the list of packages.  Hyphens or colons between the
        word ``optional`` and the first package name are allowed.  Therefore,
        you should not write ``# optional - depends on package bliss`` but simply
        ``# optional - bliss``.

      - Optional tags are case-insensitive, so you could also write ``# optional -
        Bliss``.

  If ``# optional`` or ``# needs`` is placed right after the ``sage:`` prompt,
  it is a block-scoped tag, which applies to all doctest lines until
  a blank line is encountered.

  These tags can also be applied to an entire file. If one of the first 10 lines
  of a file starts with any of ``r""" sage.doctest: optional - FEATURE``,
  ``# sage.doctest: needs FEATURE``, or ``.. sage.doctest: optional - FEATURE``
  (in ``.rst`` files), etc., then this applies to all doctests in this file.

  When a file is skipped that was explicitly given as a command line argument,
  a warning is displayed.

  .. NOTE::

       If you add such a line to a file, you are strongly encouraged
       to add a note to the module-level documentation, saying that
       the doctests in this file will be skipped unless the
       appropriate conditions are met.

- **indirect doctest:** in the docstring of a function ``A(...)``, a line
  calling ``A`` and in which the name ``A`` does not appear should have this
  flag. This prevents ``sage --coverage <file>`` from reporting the docstring as
  "not testing what it should test".

  Use it when testing special functions like ``__repr__``, ``__add__``,
  etc. Use it also when you test the function by calling ``B`` which
  internally calls ``A``:

  .. CODE-BLOCK:: rest

      This is the docstring of an ``__add__`` method. The following
      example tests it, but ``__add__`` is not written anywhere::

          sage: 1+1 # indirect doctest
          2

- **32-bit** or **64-bit:** for tests that behave differently on 32-bit or
  64-bit machines. Note that this particular flag is to be applied on the
  **output** lines, not the input lines::

      sage: hash(2^31 + 2^13)
      8193                      # 32-bit
      2147491840                # 64-bit

Per coding style (:ref:`section-coding-python`), the magic comment
should be separated by at least 2 spaces.

For multiline doctests, the comment should appear on the first
`physical line <https://docs.python.org/3/reference/lexical_analysis.html#physical-lines>`_
of the doctest (the line with the prompt ``sage:``), not on the
continuation lines (the lines with the prompt ``....:``)::

    sage: print(ZZ.random_element())        # random
    42
    sage: for _ in range(3):                # random
    ....:     print(QQ.random_element())
    1
    1/77
    -1/2

Using ``search_src`` from the Sage prompt (or ``grep``), one can
easily find the aforementioned keywords. In the case of ``todo: not
implemented``, one can use the results of such a search to direct
further development on Sage.

.. _chapter-testing:

Running automated doctests
==========================

This section describes Sage's automated testing of test files of the
following types: ``.py``, ``.pyx``, ``.sage``, ``.rst``. Briefly, use
``sage -t <file>`` to test that the examples in ``<file>`` behave
exactly as claimed. See the following subsections for more
details. See also :ref:`section-docstrings` for a discussion on how to
include examples in documentation strings and what conventions to
follow. The chapter :ref:`chapter-doctesting` contains a tutorial on
doctesting modules in the Sage library.


.. _section-testpython:

Testing .py, .pyx and .sage files
---------------------------------

Run ``sage -t <filename.py>`` to test all code examples in
``filename.py``. Similar remarks apply to ``.sage`` and ``.pyx``
files:

.. CODE-BLOCK:: shell-session

      $ sage -t [--verbose] [--optional]  [files and directories ... ]

The Sage doctesting framework is based on the standard Python doctest
module, but with many additional features (such as parallel testing,
timeouts, optional tests).  The Sage doctester recognizes ``sage:``
prompts as well as ``>>>`` prompts.  It also preparses the doctests,
just like in interactive Sage sessions.

Your file passes the tests if the code in it will run when entered
at the ``sage:`` prompt with no extra imports. Thus users are
guaranteed to be able to exactly copy code out of the examples you
write for the documentation and have them work.

For more information, see :ref:`chapter-doctesting`.


Testing reST documentation
--------------------------

Run ``sage -t <filename.rst>`` to test the examples in verbatim
environments in reST documentation.

Of course in reST files, one often inserts explanatory texts between
different verbatim environments. To link together verbatim
environments, use the ``.. link`` comment. For example:

.. CODE-BLOCK:: rest

    EXAMPLES::

            sage: a = 1


    Next we add 1 to ``a``.

    .. link::

            sage: 1 + a
            2

If you want to link all the verbatim environments together, you can
put ``.. linkall`` anywhere in the file, on a line by itself.  (For
clarity, it might be best to put it near the top of the file.)  Then
``sage -t`` will act as if there were a ``.. link`` before each
verbatim environment.  The file
:sage_root:`src/doc/en/tutorial/interfaces.rst` contains a
``.. linkall`` directive, for example.

You can also put ``.. skip`` right before a verbatim environment to
have that example skipped when testing the file.  This goes in the
same place as the ``.. link`` in the previous example.

See the files in :sage_root:`src/doc/en/tutorial/` for many
examples of how to include automated testing in reST documentation for
Sage.


.. _section-coding-general-whitespace:

General coding style regarding whitespace
=========================================

Use spaces instead of tabs for indentation. The only exception is for
makefiles, in which tabs have a syntactic meaning different from
spaces.

Do not add trailing whitespace.

Sage provides editor configuration for Emacs, using the file
``.dir-locals.el``, to use spaces instead of tabs.  Regarding trailing
whitespace, see https://www.emacswiki.org/emacs/DeletingWhitespace
for various solutions.

If you use another editor, we recommend to configure it so you do not
add tabs to files. See :ref:`section-ide`.


Global options
==============

Global options for classes can be defined in Sage using
:class:`~sage.structure.global_options.GlobalOptions`.

Miscellaneous minor things
==========================

Some decisions are arbitrary, but common conventions make life easier.

* Non-ASCII names in identifiers:

  * Translate *ä* and *ö* to *ae* and *oe*, like ``moebius_function``
    for Möbius function.
  * Translate *á* to *a*, like ``lovasz_number`` for Lovász number.

* Common function keyword arguments:

  This is a list of some keyword arguments that many functions and
  methods take.  For consistency, you should use the keywords from the
  list below with the meaning as explained here. Do not use a
  different keyword with the same meaning (for example, do not use
  ``method``; use ``algorithm`` instead).

  * ``algorithm``, a string or ``None``: choose between various
    implementation or algorithm. Use ``None`` as a default that
    selects a sensible default, which could depend on installed
    optional packages.

  * ``certificate``, a Boolean with ``False`` as default: whether the
    function should return some kind of certificate together with the
    result. With ``certificate=True`` the return value should be a
    pair `(r, c)` where `r` is the result that would be given with
    ``certificate=False`` and `c` is the certificate or ``None`` if
    there is no meaningful certificate.

  * ``proof``, a Boolean with ``True`` as default: if ``True``,
    require a mathematically proven computation. If ``False``, a
    probabilistic algorithm or an algorithm relying to non-proved
    hypothesis like RH can be used.

  * ``check``, a Boolean: do some additional checks to verify the
    input parameters. This should not otherwise influence the
    functioning of the code: if code works with ``check=True``, it should
    also work with ``check=False``.

  * ``coerce``, a Boolean: convert the input parameters to a suitable
    parent. This is typically used in constructors. You can call a
    method with ``coerce=False`` to skip some checks if the parent is
    known to be correct.

  * ``inplace``, a Boolean: whether to modify the object in-place or
    to return a copy.
