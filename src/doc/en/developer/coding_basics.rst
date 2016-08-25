.. _chapter-code-basics:

===================
General Conventions
===================


There are many ways to contribute to Sage including sharing scripts
and Sage worksheets that implement new functionality using Sage,
improving to the Sage library, or to working on the many underlying
libraries distributed with Sage [1]_.
This guide focuses on editing the Sage library itself.

Sage is not just about gathering together functionality. It is about
providing a clear, systematic and consistent way to access a large
number of algorithms, in a coherent framework that makes sense
mathematically. In the design of Sage, the semantics of objects, the
definitions, etc., are informed by how the corresponding objects are
used in everyday mathematics.

.. [1]
   See http://www.sagemath.org/links-components.html for a full list
   of packages shipped with every copy of Sage

To meet the goal of making Sage easy to read, maintain, and improve,
all Python/Cython code that is included with Sage should adhere to the
style conventions discussed in this chapter.


.. _section-coding-python:

Python Code Style
=================

Follow the standard Python formatting rules when writing code for
Sage, as explained at the following URLs:

* http://www.python.org/dev/peps/pep-0008
* http://www.python.org/dev/peps/pep-0257

In particular,

- Use 4 spaces for indentation levels. Do not use tabs as they can
  result in indentation confusion. Most editors have a feature that
  will insert 4 spaces when the tab key is hit. Also, many editors
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

      class SomeValue(object):
          def __init__(self, x):
          self._x  = 1

  and factory functions that mimic object constructors, for example
  ``PolynomialRing`` or::

       def SomeIdentityValue(x):
           return SomeValue(1)



.. _chapter-directory-structure:

Files and Directory Structure
=============================

Roughly, the Sage directory tree is layout like this. Note that we use
``SAGE_ROOT`` in the following as a shortcut for the (arbitrary) name
of the directory containing the Sage sources::

    SAGE_ROOT/
        sage          # the Sage launcher
        Makefile      # top level Makefile
        build/        # sage's build system
            deps
            install
            ...
            pkgs/     # install, patch, and metadata from spkgs
        src/
            setup.py
            module_list.py
            ...
            sage/     # sage library (formerly devel/sage-main/sage)
            ext/      # extra sage resources (formerly devel/ext-main)
            mac-app/  # would no longer have to awkwardly be in extcode
            bin/      # the scripts in local/bin that are tracked
        upstream/     # tarballs of upstream sources
        local/        # installed binaries

Python Sage library code goes into ``src/`` and uses the following
conventions. Directory names may be plural (e.g. ``rings``) and file
names are almost always singular (e.g. ``polynomial_ring.py``). Note
that the file ``polynomial_ring.py`` might still contain definitions
of several different types of polynomial rings.

.. NOTE::

   You are encouraged to include miscellaneous notes, emails, design
   discussions, etc., in your package.  Make these plain text files
   (with extension ``.txt``) in a subdirectory called ``notes``.  For
   example, see ``SAGE_ROOT/src/sage/ext/notes/``.

If you want to create a new directory in the Sage library
``SAGE_ROOT/src/sage`` (say, ``measure_theory``), that directory
should contain a file ``__init__.py`` that contains the single line
``import all`` in addition to whatever
files you want to add (say, ``borel_measure.py`` and
``banach_tarski.py``), and also a file ``all.py`` listing imports from
that directory that are important enough to be in the Sage’s global
namespace at startup.
The file ``all.py`` might look like this::

    from borel_measure import BorelMeasure
    from banach_tarski import BanachTarskiParadox

but it is generally better to use the lazy import framework::

    from sage.misc.lazy_import import lazy_import
    lazy_import('sage.measure_theory.borel_measue', 'BorelMeasure')
    lazy_import('sage.measure_theory.banach_tarski', 'BanachTarskiParadox')

Then in the file ``SAGE_ROOT/src/sage/all.py``, add a line ::

    from sage.measure_theory.all import *


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


Headings of Sage Library Code Files
===================================

The top of each Sage code file should follow this format::

    r"""
    <Very short 1-line summary>

    <Paragraph description>

    AUTHORS:

    - YOUR NAME (2005-01-03): initial version

    - person (date in ISO year-month-day format): short desc

    EXAMPLES::

    <Lots and lots of examples>
    """

    #*****************************************************************************
    #       Copyright (C) 2013 YOUR NAME <your email>
    #
    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 2 of the License, or
    # (at your option) any later version.
    #                  http://www.gnu.org/licenses/
    #*****************************************************************************

As an example, see ``SAGE_ROOT/src/sage/rings/integer.pyx`` which
contains the implementation for `\ZZ`. The ``AUTHORS:`` section is
redundant, the authoritative log for who wrote what is always the git
repository (see the output of ``git blame``). Nevertheless, it is
sometimes useful to have a very rough overview over the history,
especially if a lot of people have been working on that source file.

All code included with Sage must be licensed under the GPLv2+ or a
compatible, that is, less restrictive license (e.g. the BSD license).


.. _section-docstrings:

Documentation Strings
=====================

.. _section-docstring-function:

The docstring of a function: content
-------------------------------------

**Every** function must have a docstring that includes the following
information. You can use the existing functions of Sage as templates.

-  A **one-sentence description** of the function.

   It must be followed by a blank line and end in a period.  It
   describes the function or method's effect as a command ("Do this",
   "Return that"), not as a description like "Returns the pathname ...".

-  An **INPUT** and an **OUTPUT** block describing the input/output of
   the function. This is not optional.

   The INPUT block describes all arguments that the function accepts,
   and the OUTPUT section describes its expected output.

   1. The type names should be descriptive, but do not have to represent
      the exact Sage/Python types. For example, use "integer" for
      anything that behaves like an integer, rather than ``int``.

   2. Mention the default values of the input arguments when applicable.

   Example::

       INPUT:

       - ``p`` -- (default: 2) a positive prime integer.

       OUTPUT:

       A 5-tuple consisting of integers in this order:

       1. the smallest primitive root modulo p
       2. the smallest prime primitive root modulo p
       3. the largest primitive root modulo p
       4. the largest prime primitive root modulo p
       5. total number of prime primitive roots modulo p

   You can start the OUTPUT block with a dash if you prefer::

       OUTPUT:

       - The plaintext resulting from decrypting the ciphertext ``C``
         using the Blum-Goldwasser decryption algorithm.

-  An **EXAMPLES** block for examples. This is not optional.

   These examples are used both for:

   1. Documentation
   2. Automatic testing before each release.

   They should have good coverage of the functionality in question.

-  A **SEEALSO** block (highly recommended) with links to related parts of
   Sage. This helps users find the features that interests them and discover the
   new ones. ::

       .. SEEALSO::

           :ref:`chapter-sage_manuals_links`,
           :meth:`sage.somewhere.other_useful_method`,
           :mod:`sage.some.related.module`.

   See :ref:`chapter-sage_manuals_links` for details on how to setup
   link in Sage.

-  An **ALGORITHM** block (optional).

   It indicates what algorithm and/or what software is used, e.g.
   ``ALGORITHM: Uses Pari``. Here's a longer example with a
   bibliographical reference::

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

       REFERENCES:

       .. [Nat2000] \M. B. Nathanson. Elementary Methods in Number Theory.
          Springer, 2000.

-  A **NOTE** block for tips/tricks (optional). ::

       .. NOTE::

           You should note that this sentence is indented at least 4
           spaces. Never use the tab character.

- A **WARNING** block for critical information about your code (optional).

  For example known situations for which the code breaks, or anything
  that the user should be aware of. ::

      .. WARNING::

          Whenever you edit the Sage documentation, make sure that
          the edited version still builds. That is, you need to ensure
          that you can still build the HTML and PDF versions of the
          updated documentation. If the edited documentation fails to
          build, it is very likely that you would be requested to
          change your patch.

- A **TODO** block for future improvements (optional).

  It can contain disabled doctests to demonstrate the desired
  feature. Here's an example of a TODO block::

      .. TODO::

          Add to ``have_fresh_beers`` an interface with the faster
          algorithm "Buy a Better Fridge" (BaBF)::

              sage: have_fresh_beers('Bière de l\'Yvette', algorithm="BaBF") # not implemented
              Enjoy !

- A **PLOT** block to illustrate with pictures the output of a function.

  Generate with Sage code an object ``g`` with a ``.plot`` method, then call
  ``sphinx_plot(g)``::

      .. PLOT::

          g = graphs.PetersenGraph()
          sphinx_plot(g)

- A **REFERENCES** block to list related books or papers (optional)

  It should cite the books/research papers relevant to the code, e.g. the source
  of the algorithm that it implements. ::

      This docstring is referencing [SC]_. Just remember that references
      are global, so we can also reference to [Nat2000]_ in the ALGORITHM
      block, even if it is in a separate file. However we would not
      include the reference here since it would cause a conflict.

      REFERENCES:

      .. [SC] Conventions for coding in sage.
         http://doc.sagemath.org/html/en/developer/conventions.html.

  When abbreviating the first name of an author, be sure to put a
  backslash in front of it. This ensures that the letter (``C.``
  in the example below) will not be interpreted as a list
  enumerator::

      .. [Gauss] \C. F. Gauss, Disquisitiones Arithmeticae, 1801.

  See the `Sphinx/ReST markup for citations <http://sphinx.pocoo.org/rest.html#citations>`_. For links toward trac tickets or wikipedia, see :ref:`chapter-sage_manuals_links`.

- A **TESTS** block (optional)

  Formatted just like EXAMPLES, containing tests that are not relevant
  to users.  In particular, these blocks are not shown when users ask
  for help via `foo?`: they are stripped by the function
  :func:`sage.misc.sagedoc.skip_TESTS_block`.

  For the purposes of removal, A "TESTS" block is a block starting
  with "TEST:" or "TESTS:" (or the same with two colons), on a line on
  its own, and ending with an unindented line (that is, the same level
  of indentation as "TESTS") matching one of the following:

  one of the following ways:

  - a line which starts with whitespace and then a Sphinx directive
    of the form ".. foo:", optionally followed by other text.

  - a line which starts with whitespace and then text of the form
    "UPPERCASE:", optionally followed by other text.

  - lines which look like a ReST header: one line containing
    anything, followed by a line consisting only of whitespace,
    followed by a string of hyphens, equal signs, or other
    characters which are valid markers for ReST headers:
    ``- = ` : ' " ~ _ ^ * + # < >``.

Template
^^^^^^^^

Use the following template when documenting functions. Note the
indentation:

.. skip    # do not doctest

::

    def point(self, x=1, y=2):
        r"""
        Return the point `(x^5,y)`.

        INPUT:

        - ``x`` -- integer (default: 1) the description of the
          argument ``x`` goes here.  If it contains multiple lines, all
          the lines after the first need to begin at the same indentation
          as the backtick.

        - ``y`` -- integer (default: 2) the ...

        OUTPUT:

        The point as a tuple.

        .. SEEALSO::

            :func:`line`

        EXAMPLES:

        This example illustrates ...

        ::

            sage: A = ModuliSpace()
            sage: A.point(2,3)
            xxx

        We now ...

        ::

            sage: B = A.point(5,6)
            sage: xxx

        It is an error to ...::

            sage: C = A.point('x',7)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'r' to an integer

        .. NOTE::

            This function uses the algorithm of [BCDT]_ to determine
            whether an elliptic curve `E` over `Q` is modular.

        ...

        REFERENCES:

        .. [BCDT] Breuil, Conrad, Diamond, Taylor,
           "Modularity ...."
        """
        <body of the function>

You are strongly encouraged to:

- Use LaTeX typesetting (see :ref:`section-latex-typeset`).

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
there is not one already. That is, you can do the following::

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

LaTeX Typesetting
-----------------

In Sage's documentation LaTeX code is allowed and is marked with **backticks or
dollar signs**:

    ```x^2 + y^2 = 1``` and ``$x^2 + y^2 = 1$`` both yield `x^2 + y^2 = 1`.

**Backslashes:** For LaTeX commands containing backslashes, either use double
backslashes or begin the docstring with a ``r"""`` instead of ``"""``. Both of
the following are valid::

    def cos(x):
        """
        Return `\\cos(x)`.
        """

    def sin(x):
        r"""
        Return $\sin(x)$.
        """

**MATH block:** This is similar to the LaTeX syntax ``\[<math expression>\]`` (or
``$$<math expression>$$``). For instance::

    .. MATH::

        \sum_{i=1}^{\infty} (a_1 a_2 \cdots a_i)^{1/i}
        \leq
        e \sum_{i=1}^{\infty} a_i

.. MATH::

    \sum_{i=1}^{\infty} (a_1 a_2 \cdots a_i)^{1/i}
    \leq
    e \sum_{i=1}^{\infty} a_i

The **aligned** environment works as it does in LaTeX::

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
<http://sphinx-doc.org/latest/ext/math.html?highlight=nowrap#directive-math>`_. ::

    .. MATH::
       :nowrap:

       \begin{align}
          1+...+n &= n(n+1)/2\\
          &= O(n^2)\\
       \end{tabular}

.. MATH::
   :nowrap:

   \begin{align}
   1+...+n &= n(n+1)/2\\
   &= O(n^2)\\
   \end{align}

**Readability balance:** in the interactive console, LaTeX formulas contained in
the documentation are represented by their LaTeX code (with
backslashes stripped). In this situation ``\\frac{a}{b}`` is less readable than ``a/b``
or ``a b^{-1}`` (some users may not even know LaTeX code). Make it pleasant for
everybody as much as you can manage.

**Commons rings** `(\Bold{Z},\Bold{N},...)`: The Sage LaTeX style is to typeset
standard rings and fields using the locally-defined macro ``\\Bold``
(e.g. ``\\Bold{Z}`` gives `\Bold{Z}`).

**Shortcuts** are available which preserve readability, e.g. ``\\ZZ`` (`\ZZ`),
``\\RR`` (`\RR`), ``\\CC`` (`\CC`), and ``\\QQ`` (`\QQ`). They appear as
LaTeX-formatted ``\\Bold{Z}`` in the html manual, and as ``Z`` in the
interactive help. Other examples: ``\\GF{q}`` (`\GF{q}`) and ``\\Zmod{p}``
(`\Zmod{p}`).

See the file ``SAGE_ROOT/src/sage/misc/latex_macros.py`` for a full list and
for details about how to add more macros.

.. _section-doctest-writing:

Writing Testable Examples
-------------------------

The examples from Sage's documentation have a double purpose:

- They provide **illustrations** of the code's usage to the users

- They are **tests** that are checked before each release, helping us avoid
  new bugs.

All new doctests added to Sage should **pass all tests** (see
:ref:`chapter-doctesting`), i.e. running ``sage -t your_file.py`` should not
give any error messages. Below are instructions about how doctests should be
written.

**What doctests should test:**

- **Interesting examples** of what the function can do. This will be the
  most helpful to a lost user. It is also the occasion to check famous
  theorems (just in case)::

    sage: is_prime(6) # 6 is not prime
    False
    sage: 2 * 3 # and here is a proof
    6

- All **meaningful combinations** of input arguments. For example a function may
  accept an ``algorithm="B"`` argument, and doctests should involve both
  ``algorithm="A"`` and ``algorithm="B"``.

- **Corner cases:** the code should be able to handle a 0 input, or an empty
  set, or a null matrix, or a null function, ... All corner cases should be
  checked, as they are the most likely to be broken, now or in the future. This
  probably belongs to the TESTS block (see :ref:`section-docstring-function`).

- **Systematic tests** of all small-sized inputs, or tests of **random**
  instances if possible.

  .. NOTE::

     Note that **TestSuites** are an automatic way to generate some of these
     tests in specific situations. See
     ``SAGE_ROOT/src/sage/misc/sage_unittest.py``.

**The syntax:**

- **Environment:** doctests should work if you copy/paste them in Sage's
  interactive console. For example, the function ``AA()`` in the file
  ``SAGE_ROOT/src/sage/algebras/steenrod/steenrod_algebra.py`` includes an
  EXAMPLES block containing the following::

    sage: from sage.algebras.steenrod.steenrod_algebra import AA as A
    sage: A()
    mod 2 Steenrod algebra, milnor basis

  Sage does not know about the function ``AA()`` by default, so it needs to be
  imported before it is tested. Hence the first line in the example.

- **Preparsing:** As in Sage's console, `4/3` returns `4/3` and not `1` as in
  Python 2.7. Testing occurs with full Sage preparsing of input within the standard
  Sage shell environment, as described in :ref:`section-preparsing`.

- **Writing files:** If a test outputs to a file, the file should be a temporary
  file.  Use :func:`tmp_filename` to get a temporary filename, or
  :func:`tmp_dir` to get a temporary directory. An example from
  ``SAGE_ROOT/src/sage/plot/graphics.py``)::

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

Special Markup to Influence Tests
---------------------------------

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

  However, most functions generating pseudorandom output do not need this tag
  since the doctesting framework guarantees the state of the pseudorandom number
  generators (PRNGs) used in Sage for a given doctest.

  When possible, avoid the problem, e.g.: rather than checking the value of the
  hash in a doctest, one could illustrate successfully using it as a key in a
  dict.

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
      sage: e^(i*pi/4).n() # rel tol 1e-1
      0.7 + 0.7*I
      sage: ((x+1.001)^4).expand() # rel tol 2
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

      sage: factor(x*y - x*z)    # todo: not implemented

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

- **optional:** A line flagged with ``optional - keyword`` is not tested unless
  the ``--optional=keyword`` flag is passed to ``sage -t`` (see
  :ref:`section-optional-doctest-flag`). The main applications are:

  - **optional packages:** When a line requires an optional package to be
    installed (e.g. the ``sloane_database`` package)::

      sage: SloaneEncyclopedia[60843]    # optional - sloane_database

  - **internet:** For lines that require an internet connection::

       sage: sloane_sequence(60843)       # optional - internet

  - **bug:** For lines that describe bugs. Alternatively, use ``# known bug``
    instead: it is an alias for ``optional bug``. ::

        The following should yield 4.  See :trac:`2`. ::

            sage: 2+2  # optional: bug
            5
            sage: 2+2  # known bug
            5

  .. NOTE::

      - Any words after ``# optional`` are interpreted as a list of
        package names, separated by spaces.

      - Any punctuation (periods, commas, hyphens, semicolons, ...) after the
        first word ends the list of packages.  Hyphens or colons between the
        word ``optional`` and the first package name are allowed.  Therefore,
        you should not write ``optional: needs package CHomP`` but simply
        ``optional: CHomP``.

      - Optional tags are case-insensitive, so you could also write ``optional:
        chOMP``.

- **indirect doctest:** in the docstring of a function ``A(...)``, a line
  calling ``A`` and in which the name ``A`` does not appear should have this
  flag. This prevents ``sage --coverage <file>`` from reporting the docstring as
  "not testing what it should test".

  Use it when testing special functions like ``__repr__``, ``__add__``,
  etc. Use it also when you test the function by calling ``B`` which
  internally calls ``A``::

      This is the docstring of an ``__add__`` method. The following
      example tests it, but ``__add__`` is not written anywhere::

          sage: 1+1 # indirect doctest
          2

- **32-bit** or **64-bit:** for tests that behave differently on 32-bit or
  64-bit machines. Note that this particular flag is to be applied on the
  **output** lines, not the input lines::

      sage: hash(-920390823904823094890238490238484)
      -873977844            # 32-bit
      6874330978542788722   # 64-bit

Using ``search_src`` from the Sage prompt (or ``grep``), one can
easily find the aforementioned keywords. In the case of ``todo: not
implemented``, one can use the results of such a search to direct
further development on Sage.

.. _chapter-testing:

Running Automated Tests
=======================

This section describes Sage's automated testing of test files of the
following types: ``.py``, ``.pyx``, ``.sage``, ``.rst``. Briefly, use
``sage -t <file>`` to test that the examples in ``<file>`` behave
exactly as claimed. See the following subsections for more
details. See also :ref:`section-docstrings` for a discussion on how to
include examples in documentation strings and what conventions to
follow. The chapter :ref:`chapter-doctesting` contains a tutorial on
doctesting modules in the Sage library.


.. _section-testpython:

Testing .py, .pyx and .sage Files
---------------------------------

Run ``sage -t <filename.py>`` to test all code examples in
``filename.py``. Similar remarks apply to ``.sage`` and ``.pyx``
files::

      sage -t [--verbose] [--optional]  [files and directories ... ]

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


Testing ReST Documentation
--------------------------

Run ``sage -t <filename.rst>`` to test the examples in verbatim
environments in ReST documentation.

Of course in ReST files, one often inserts explanatory texts between
different verbatim environments. To link together verbatim
environments, use the ``.. link`` comment. For example::

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
``SAGE_ROOT/src/doc/en/tutorial/interfaces.rst`` contains a
``.. linkall`` directive, for example.

You can also put ``.. skip`` right before a verbatim environment to
have that example skipped when testing the file.  This goes in the
same place as the ``.. link`` in the previous example.

See the files in ``SAGE_ROOT/src/doc/en/tutorial/`` for many
examples of how to include automated testing in ReST documentation for
Sage.

.. _chapter-picklejar:

The Pickle Jar
==============

Sage maintains a pickle jar at
``SAGE_ROOT/src/ext/pickle_jar/pickle_jar.tar.bz2`` which is a tar
file of "standard" pickles created by ``sage``. This pickle jar is
used to ensure that sage maintains backward compatibility by
having :func:`sage.structure.sage_object.unpickle_all` check that
``sage`` can always unpickle all of the pickles in the pickle jar as
part of the standard doc testing framework.

Most people first become aware of the pickle_jar when their patch breaks the
unpickling of one of the "standard" pickles in the pickle jar due to the
failure of the doctest::

    sage -t src/sage/structure/sage_object.pyx

When this happens an error message is printed which contains the following
hints for fixing the uneatable pickle::

    ----------------------------------------------------------------------
    ** This error is probably due to an old pickle failing to unpickle.
    ** See sage.structure.sage_object.register_unpickle_override for
    ** how to override the default unpickling methods for (old) pickles.
    ** NOTE: pickles should never be removed from the pickle_jar!
    ----------------------------------------------------------------------

For more details about how to fix unpickling errors in the pickle jar
see :func:`sage.structure.sage_object.register_unpickle_override`

.. WARNING::

    Sage's pickle jar helps to ensure backward compatibility in sage. Pickles should
    **only** be removed from the pickle jar after the corresponding objects
    have been properly deprecated. Any proposal to remove pickles from the
    pickle jar should first be discussed on sage-devel.


Global Options
==============

Global options for classes can be defined in Sage using
:class:`~sage.structure.global_options.GlobalOptions`.

