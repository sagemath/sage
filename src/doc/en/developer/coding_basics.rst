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


An Example Is Worth a Thousand Words
====================================

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
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.
    #                  http://www.gnu.org/licenses/
    #*****************************************************************************

As an example, see ``SAGE_ROOT/src/sage/rings/integer.pyx`` which
contains the implementation for `\ZZ`. The ``AUTHORS:`` section is
redundant, the authoritative log for who wrote what is always the git
repository (see the output of ``git blame``). Nevertheless, it is
sometimes useful to have a very rough overview over the history,
especially if a lot of people have been working on that source file.

All code included with Sage must be licensed under the GPLv3+ or a
compatible, that is, less restrictive license (e.g. GPLv2+ or the BSD
license).


.. _section-docstrings:

Documentation Strings
=====================


Docstring Markup With ReST and Sphinx
-------------------------------------

**Every** function must have a docstring that includes the following
information. Source files in the Sage library contain numerous
examples on how to format your documentation, so you could use them as
a guide.

-  A one-sentence description of the function, followed by a blank line
   and ending in a period. It prescribes the function or method's
   effect as a command ("Do this", "Return that"), not as a
   description; e.g. don't write "Returns the pathname ...".

-  An INPUT and an OUTPUT block for input and output arguments (see
   below for format). The type names should be descriptive, but do not
   have to represent the exact Sage/Python types. For example, use
   "integer" for anything that behaves like an integer; you do not have
   to put a precise type name such as ``int``. The INPUT block
   describes the expected input to your function or method, while the
   OUTPUT block describes the expected output of the
   function/method. If appropriate, you need to describe any default
   values for the input arguments. For example::

       INPUT:

       - ``p`` -- (default: 2) a positive prime integer.

       OUTPUT:

       A 5-tuple consisting of integers in this order:

       1. the smallest primitive root modulo p
       2. the smallest prime primitive root modulo p
       3. the largest primitive root modulo p
       4. the largest prime primitive root modulo p
       5. total number of prime primitive roots modulo p

   Some people prefer to format their OUTPUT section as a block by
   using a dash. That is acceptable as well::

       OUTPUT:

       - The plaintext resulting from decrypting the ciphertext ``C``
         using the Blum-Goldwasser decryption algorithm.

-  An EXAMPLES block for examples. This is not optional. These
   examples are used both for documentation and for automatic testing
   before each release so should have good coverage of the functionality
   in question. New functions without these doctests will not be accepted
   for inclusion with Sage.

-  A SEEALSO block (optional) with links to related things in Sage. A SEEALSO
   block should start with ``.. SEEALSO::``. It can also be the lower-case form
   ``.. seealso::``. However, you are encouraged to use the upper-case form
   ``.. SEEALSO::``. See :ref:`chapter-sage_manuals_links` for details on how
   to setup link in Sage.  Here's an example of a SEEALSO block::

       .. SEEALSO::

           :ref:`chapter-sage_manuals_links`

-  An ALGORITHM block (optional) which indicates what algorithm
   and/or what software is used. For example
   ``ALGORITHM: Uses Pari``. Here's a longer example that describes an
   algorithm used. Note that it also cites the reference where this
   algorithm can be found::

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

       .. [Nat2000] M.B. Nathanson. Elementary Methods in Number Theory.
          Springer, 2000.

   You can also number the steps in your algorithm using the hash-dot
   symbol. This way, the actual numbering of the steps are
   automatically taken care of when you build the documentation::

        ALGORITHM:

        The Blum-Goldwasser decryption algorithm is described in Algorithm
        8.56, page 309 of [MenezesEtAl1996]_. The algorithm works as follows:

        #. Let `C` be the ciphertext `C = (c_1, c_2, \dots, c_t, x_{t+1})`.
           Then `t` is the number of ciphertext sub-blocks and `h` is the
           length of each binary string sub-block `c_i`.
        #. Let `(p, q, a, b)` be the private key whose corresponding
           public key is `n = pq`. Note that `\gcd(p, q) = ap + bq = 1`.
        #. Compute `d_1 = ((p + 1) / 4)^{t+1} \bmod{(p - 1)}`.
        #. Compute `d_2 = ((q + 1) / 4)^{t+1} \bmod{(q - 1)}`.
        #. Let `u = x_{t+1}^{d_1} \bmod p`.
        #. Let `v = x_{t+1}^{d_2} \bmod q`.
        #. Compute `x_0 = vap + ubq \bmod n`.
        #. For `i` from 1 to `t`, do:

           #. Compute `x_i = x_{t-1}^2 \bmod n`.
           #. Let `p_i` be the `h` least significant bits of `x_i`.
           #. Compute `m_i = p_i \oplus c_i`.

        #. The plaintext is `m = m_1 m_2 \cdots m_t`.

-  A NOTE block for special notes (optional). Include information
   such as purpose etc. A NOTE block should start with
   ``.. NOTE::``. You can also use the lower-case version
   ``.. note::``, but do not mix lower-case with upper-case. However,
   you are encouraged to use the upper-case version ``.. NOTE::``. If
   you want to put anything within the NOTES block, you should
   indent it at least 4 spaces (no tabs). Here's an example of a NOTE
   block::

       .. NOTE::

           You should note that this sentence is indented at least 4
           spaces. Avoid tab characters as much as possible when
           writing code or editing the Sage documentation. You should
           follow Python conventions by using spaces only.

- A WARNING block for critical information about your code. For
  example, the WARNING block might include information about when or
  under which conditions your code might break, or information that
  the user should be particularly aware of. A WARNING block should start
  with ``.. WARNING::``. It can also be the lower-case form
  ``.. warning::``. However, you are encouraged to use the upper-case
  form ``.. WARNING::``. Here's an example of a WARNING block::

      .. WARNING::

          Whenever you edit the Sage documentation, make sure that
          the edited version still builds. That is, you need to ensure
          that you can still build the HTML and PDF versions of the
          updated documentation. If the edited documentation fails to
          build, it is very likely that you would be requested to
          change your patch.

- A TODO block for room for improvements. The TODO block might
  contains disabled doctests to demonstrate the desired feature.  A TODO block
  should start with ``.. TODO::``. It can also be the lower-case form
  ``.. todo::``. However, you are encouraged to use the upper-case form
  ``.. TODO::``. Here's an example of a TODO block::

      .. TODO::

          Improve further function ``have_fresh_beers`` using algorithm
          ``buy_a_better_fridge``::

              sage: have_fresh_beers('Bière de l\'Yvette') # todo: not implemented
              Enjoy !

- A REFERENCES block to list books or papers (optional). This block serves
  a similar purpose to a list of references in a research paper, or a
  bibliography in a monograph. If your method, function or class uses an
  algorithm that can be found in a standard reference, you should list
  that reference under this block. The Sphinx/ReST markup for
  citations is described at
  http://sphinx.pocoo.org/rest.html#citations. See below for an example.
  Sage also add specific markup for links to sage trac tickets and
  Wikipedia. See :ref:`chapter-sage_manuals_links`. Here's an example of a
  REFERENCES block::

      This docstring is referencing [SC]_. Just remember that references
      are global, so we can also reference to [Nat2000]_ in the ALGORITHM
      block, even if it is in a separate file. However we would not
      include the reference here since it would cause a conflict.

      REFERENCES:

      .. [SC] Conventions for coding in sage.
         http://www.sagemath.org/doc/developer/conventions.html.

- A TESTS block (optional), formatted just like EXAMPLES, for additional
  tests which should be part of the regression suite but are not
  illustrative enough to merit placement in EXAMPLES.

Use the following template when documenting functions. Note the
indentation

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
            TypeError: unable to convert x (=r) to an integer

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

- Use nice LaTeX formatting everywhere, see
  :ref:`section-latex-typeset`.

- Liberally describe what the examples do. Note that there must be a
  blank line after the example code and before the explanatory text
  for the next example (indentation is not enough).

- Illustrate any exceptions raised by the function with examples, as
  given above. (It is an error to ...; In particular, use ...)

- Include many examples. These are automatically tested on a regular
  basis, and are crucial for the quality and adaptability of
  Sage. Without such examples, small changes to one part of Sage that
  break something else might not go seen until much later when someone
  uses the system, which is unacceptable. Note that new functions
  without doctests will not be accepted for inclusion in Sage.

Functions whose names start with an underscore are considered
private. Hence they do not appear in the reference manual, and their
docstring should not contain any information that is crucial for Sage
users. Having said that, you can explicitly enable their docstrings to
be shown as part of the documentation of another method. For example::

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

An EXAMPLES or TESTS block is still required for these private functions.

A special case is the constructor ``__init__``, which clearly starts
with an underscore. However, due to its special status the
``__init__`` docstring is used as the class docstring if there is not
one already. That is, you can do the following::

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


In ReST documentation, you use backticks \` to mark LaTeX code to be
typeset.  In Sage docstrings, you may also use dollar signs instead.
Thus ```x^2 + y^2 = 1``` and ``$x^2 + y^2 = 1$`` should produce
identical output. If you use TeX commands containing backslashes in
docstrings, then either use double backslashes or place an "r" right
before the first triple opening quote. For example, both of the
following are valid::

    def cos(x):
        """
        Return `\\cos(x)`.
        """

    def sin(x):
        r"""
        Return $\sin(x)$.
        """

You can also use the MATH block to format complicated mathematical
expressions::

    .. MATH::

        \sum_{i=1}^{\infty} (a_1 a_2 \cdots a_i)^{1/i}
        \leq
        e \sum_{i=1}^{\infty} a_i

Note that the MATH block is automatically wrapped in a latex math
environment (i.e. in ``\[ \]`` or ``$$``, etc.). To use aligned equations,
use the **aligned** environment::

    .. MATH::

        \begin{aligned}
         f(x) & = x^2 - 1 \\
         g(x) & = x^x - f(x - 2)
        \end{aligned}

If you wish to explicitly not wrap the MATH block, make the first line of
the indented block ``:nowrap:``::

    .. MATH::
        :nowrap:

        This is now plain text so I can do things like $x = 5$.

.. WARNING::

    With or without ``:nowrap:``, the *html* documentation output
    currently will work if you use environments such as **align**
    which wrap their contents in math mode. However, ``:nowrap:``
    is necessary for the *pdf* documentation to build correctly.

The Sage LaTeX style is to typeset standard rings and fields like the
integers and the real numbers using the locally-defined macro
``\\Bold``, as in ``\\Bold{Z}`` for the integers. This macro is
defined to be ordinary bold-face ``\\mathbf`` by default, but users
can switch to blackboard-bold ``\\mathbb`` and back on-the-fly by
using ``latex.blackboard_bold(True)`` and
``latex.blackboard_bold(False)``.

The docstring will be available interactively (for the "def point..."
example above, by typing "point?" at the "sage:" prompt) and also in
the reference manual. When viewed interactively, LaTeX code has the
backslashes stripped from it, so "\\cos" will appear as "cos".

Because of the dual role of the docstring, you need to strike a
balance between readability (for interactive help) and using perfect
LaTeX code (for the reference manual).  For instance, instead of using
"\\frac{a}{b}", use "a/b" or maybe "a b^{-1}".  Also keep in mind that
some users of Sage are not familiar with LaTeX; this is another reason
to avoid complicated LaTeX expressions in docstrings, if at all
possible: "\\frac{a}{b}" will be obscure to someone who doesn't know
any LaTeX.

Finally, a few non-standard LaTeX macros are available to help achieve
this balance, including "\\ZZ", "\\RR", "\\CC", and "\\QQ".  These are
names of Sage rings, and they are typeset using a single boldface
character; they allow the use of "\\ZZ" in a docstring, for example,
which will appear interactively as "ZZ" while being typeset as
"\\Bold{Z}" in the reference manual.  Other examples are "\\GF" and
"\\Zmod", each of which takes an argument: "\\GF{q}" is typeset as
"\\Bold{F}_{q}" and "\\Zmod{n}" is typeset as "\\Bold{Z}/n\\Bold{Z}".
See the file ``$SAGE_ROOT/src/sage/misc/latex_macros.py`` for a
full list and for details about how to add more macros.


Writing Testable Examples
-------------------------

The code in the examples should pass automatic testing. This means
that if the above code is in the file ``f.py`` (or ``f.sage``), then
``sage -t f.py`` should not give any error messages. Testing occurs
with full Sage preparsing of input within the standard Sage shell
environment, as described in :ref:`section-preparsing`. **Important:**
The file ``f.py`` is not imported when running tests unless you have
arranged that it be imported into your Sage environment, i.e. unless
its functions are available when you start Sage using the ``sage``
command. For example, the function ``AA()`` in the file
``SAGE_ROOT/src/sage/algebras/steenrod/steenrod_algebra.py`` includes
an EXAMPLES block containing the following::

    sage: from sage.algebras.steenrod.steenrod_algebra import AA as A
    sage: A()
    mod 2 Steenrod algebra, milnor basis

Sage does not know about the function ``AA()`` by default, so
it needs to be imported before it is tested. Hence the first line in
the example.

When writing documentation, keep the following points in mind:

- All input is preparsed before being passed to Python, e.g. ``2/3``
  is replaced by ``Integer(2)/Integer(3)``, which evaluates to ``2/3``
  as a rational instead of the Python int ``0``. For more information
  on preparsing, see :ref:`section-preparsing`.

- If a test outputs to a file, the file should be a temporary file.
  Use :func:`tmp_filename` to get a temporary filename, or
  :func:`tmp_dir` to get a temporary directory.  For example (taken
  from the file ``SAGE_ROOT/src/sage/plot/graphics.py``)::

      sage: plot(x^2 - 5, (x, 0, 5), ymin=0).save(tmp_filename(ext='.png'))

- You may write tests that span multiple lines.  The best way to do so
  is to use the line continuation marker ``....:`` ::

      sage: for n in srange(1,10):
      ....:     if n.is_prime():
      ....:         print n,
      2 3 5 7

  If you have a long line of code, you may want to consider adding a
  backslash to the end of the line, which tells the doctesting
  framework to join that current line with the next.  This syntax is
  non-standard so may be removed in a future version of Sage, but in
  the mean time it can be useful for breaking up large integers across
  multiple lines::

      sage: n = 123456789123456789123456789\
      ....:     123456789123456789123456789
      sage: n.is_prime()
      False


.. _section-further_conventions:

Special Markup to Influence Tests
---------------------------------

There are a number of magic comments that you can put into the example
code that change how the output is verified by the Sage doctest
framework. Here is a comprehensive list:

- If a test line contains the comment ``random``, it is executed but
  it is not checked that the output agrees with the output in the
  documentation string. For example, the docstring for the
  ``__hash__`` method for ``CombinatorialObject`` in
  ``SAGE_ROOT/src/sage/combinat/combinat.py`` includes the lines::

      sage: c = CombinatorialObject([1,2,3])
      sage: hash(c)   # random
      1335416675971793195
      sage: c.__hash__()   # random
      1335416675971793195

  However, most functions generating pseudorandom output do not need
  this tag since the doctesting framework guarantees the state of the
  pseudorandom number generators (PRNGs) used in Sage for a given
  doctest. See :ref:`chapter-randomtesting` for details on this
  framework.  It is preferable to write tests that do not expose this
  non-determinism, for example rather than checking the value of the
  hash in a dockets, one could illustrate successfully using it as a
  key in a dict.

- If a line contains the comment ``long time`` then that line is not
  tested unless the ``--long`` option is given, e.g.  ``sage -t --long
  f.py``. Use this to include examples that take more than about a
  second to run. These will not be run regularly during Sage
  development, but will get run before major releases. No example
  should take more than about 30 seconds.

  For instance, here is part of the docstring from the ``regulator``
  method for rational elliptic curves, from the file
  ``SAGE_ROOT/devel/sage/sage/schemes/elliptic_curves/ell_rational.py``::

      sage: E = EllipticCurve([0, 0, 1, -1, 0])
      sage: E.regulator()        # long time (1 second)
      0.0511114082399688

- If a comment contains ``tol`` or ``tolerance``, numerical results are
  only verified to the given tolerance. This may be prefixed by
  ``abs[olute]`` or ``rel[ative]`` to specify whether to measure
  absolute or relative error; this defaults to relative error except
  when the expected value is exactly zero::

      sage: RDF(pi)                               # abs tol 1e-5
      3.14159
      sage: [10^n for n in [0.0 .. 4]]            # rel tol 2e-4
      [0.9999, 10.001, 100.01, 999.9, 10001]

  This can be useful when the exact output is subject to rounding
  error and/or processor floating point arithmetic variation.  Here
  are some more examples.

  A singular value decomposition of a matrix will produce two unitary
  matrices.  Over the reals, this means the inverse of the matrix is
  equal to its transpose.  We test this result by applying the norm to
  a matrix difference.  The result will usually be a "small" number,
  distinct from zero::

      sage: A = matrix(RDF, 8, range(64))
      sage: U, S, V = A.SVD()
      sage: (U.transpose()*U-identity_matrix(8)).norm(p=2)    # abs tol 1e-10
      0.0

  The 8-th cyclotomic field is generated by the complex number
  `e^\frac{i\pi}{4}`.  Here we compute a numerical approximation::

      sage: K.<zeta8> = CyclotomicField(8)
      sage: N(zeta8)                             # absolute tolerance 1e-10
      0.7071067812 + 0.7071067812*I

  A relative tolerance on a root of a polynomial.  Notice that the
  root should normally print as ``1e+16``, or something similar.
  However, the tolerance testing causes the doctest framework to use
  the output in a *computation*, so other valid text representations
  of the predicted value may be used.  However, they must fit the
  pattern defined by the regular expression ``float_regex`` in
  :mod:`sage.doctest.parsing`::

      sage: y = polygen(RDF, 'y')
      sage: p = (y - 10^16)*(y-10^(-13))*(y-2); p
      y^3 - 1e+16*y^2 + 2e+16*y - 2000.0
      sage: p.roots(multiplicities=False)[2]     # relative tol 1e-10
      10000000000000000

- If a comment contains ``not implemented`` or ``not tested``, it is
  never tested. It is good to include lines like this to make clear
  what we want Sage to eventually implement::

      sage: factor(x*y - x*z)    # todo: not implemented

  It is also immediately clear to the user that the indicated example
  does not currently work.

- If one of the first 10 lines of a file starts with ``r"""
  nodoctest`` (or ``""" nodoctest`` or ``# nodoctest`` or ``%
  nodoctest`` or ``.. nodoctest``, or any of these with different
  spacing), then that file will be skipped.  If a directory contains a
  file ``nodoctest.py``, then that whole directory will be
  skipped. Neither of this applies to files or directories which are
  explicitly given as command line arguments: those are always tested.

- If a comment contains ``optional - PKGNAME``, it is not tested
  unless the ``--optional=PKGNAME`` flag is passed to ``sage -t``.
  Mark a doctest as ``optional`` if it requires optional packages.
  Running ``sage -t --optional=all f.py`` executes all doctests,
  including all optional tests.  Running
  ``sage -t --optional=sage,sloane_database f.py`` runs the normal
  tests (because of ``--optional=sage``), as well as those marked as
  ``# optional - sloane_database``.  For example, the file
  ``SAGE_ROOT/src/sage/databases/sloane.py`` contains the lines::

       sage: sloane_sequence(60843)       # optional - internet

  and::

       sage: SloaneEncyclopedia[60843]    # optional - sloane_database

  The first of these just needs internet access, while the second
  requires that the "sloane_database" package be installed.  Calling
  ``sage -t --optional=all`` on this file runs both of these tests,
  while calling ``sage -t --optional=sage,internet`` on it will only
  run the first test.  A test requiring several packages should be
  marked ``# optional - pkg1 pkg2`` and executed by
  ``sage -t --optional=sage,pkg1,pkg2 f.py``.

  .. NOTE::

      Any words after ``# optional`` are interpreted as a list of
      package names, separated by spaces.  Any punctuation (periods,
      commas, hyphens, semicolons, ...)  after the first word ends the
      list of packages.  Hyphens or colons between the word
      ``optional`` and the first package name are allowed.  Therefore,
      you should not write ``optional: needs package CHomP`` but
      simply ``optional: CHomP``.  Optional tags are case-insensitive,
      so you could also write ``optional: cHoMp``.

- If you are documenting a known bug in Sage, mark it as ``known bug``
  or ``optional: bug``.  For example::

      The following should yield 4.  See :trac:`2`. ::

          sage: 2+2  # optional: bug
          5

  Then the doctest will be skipped by default, but could be revealed
  by running ``sage -t --optional=sage,bug ...``.  (A doctest marked
  as ``known bug`` gets automatically converted to ``optional bug``).

- Some tests (hashing for example) behave differently on 32-bit and
  64-bit platforms.  You can mark a line (generally the output) with
  either ``# 32-bit`` or ``# 64-bit`` and the testing framework will
  remove any lines that don't match the current architecture.  For
  example::

      sage: z = 32
      sage: z.powermodm_ui(2^32-1, 14)
      Traceback (most recent call last):                              # 32-bit
      ...                                                             # 32-bit
      OverflowError: exp (=4294967295) must be <= 4294967294          # 32-bit
      8              # 64-bit

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
``SAGE_ROOT/devel/sage/doc/en/tutorial/interfaces.rst`` contains a
``.. linkall`` directive, for example.

You can also put ``.. skip`` right before a verbatim environment to
have that example skipped when testing the file.  This goes in the
same place as the ``.. link`` in the previous example.

See the files in ``SAGE_ROOT/devel/sage/doc/en/tutorial/`` for many
examples of how to include automated testing in ReST documentation for
Sage.

.. _chapter-picklejar:

The Pickle Jar
==============

Sage maintains a pickle jar at
``SAGE_ROOT/src/ext/pickle_jar/pickle_jar.tar.bz2`` which is a tar
file of "standard" pickles created by ``sage``. This pickle jar is
used to ensure that sage maintains backward compatibility by have
having :func:`sage.structure.sage_object.unpickle_all` check that
``sage`` can always unpickle all of the pickles in the pickle jar as
part of the standard doc testing framework.

Most people first become aware of the pickle_jar when their patch breaks the
unpickling of one of the "standard" pickles in the pickle jar due to the
failure of the doctest::

    sage -t devel/sage-main/sage/structure/sage_object.pyx

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


.. _chapter-randomtesting:

Randomized Testing
==================

In addition to all the examples in your docstrings, which serve as
both demonstrations and tests of your code, you should consider
creating a test suite. Think of this as a program that will run for a
while and "tries" to crash your code using randomly generated
input. Your test code should define a class ``Test`` with a
``random()`` method that runs random tests. These are all assembled
together later, and each test is run for a certain amount of time on a
regular basis.

For an example, see the file
``SAGE_ROOT/src/sage/modular/modsym/tests.py``.


Global Options
==============

Global options for classes can be defined in Sage using
:class:`~sage.structure.global_options.GlobalOptions`.

