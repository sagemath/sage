.. highlight:: bash

.. _chapter-sage_manuals:

================
The Sage Manuals
================

Sage's manuals are written in `ReST <http://docutils.sourceforge.net/rst.html>`_
(reStructuredText), and generated with the software `Sphinx
<https://www.sphinx-doc.org/>`_:

.. LIST-TABLE::
   :widths: 4 12
   :header-rows: 1

   * - Name
     - Files

   * - `Tutorial <../tutorial/index.html>`_
     - :sage_root:`src/doc/en/tutorial`

   * - `Developer's guide <../developer/index.html>`_
     - :sage_root:`src/doc/en/developer`

   * - `Constructions <../constructions/index.html>`_
     - :sage_root:`src/doc/en/constructions`

   * - `Installation guide <../installation/index.html>`_
     - :sage_root:`src/doc/en/installation`

   * - `Reference manual <../reference/index.html>`_
     - :sage_root:`src/doc/en/reference`
       (most of it is generated from the
       source code)

- Additionally, more specialized manuals can be found under
  :sage_root:`src/doc/en`.

- Some documents have been **translated** into other languages. In order to
  access them, change ``en/`` into ``fr/``, ``es/``, ``de/``... See :ref:`section-manuals-names`.

.. _section-manuals-edit:

Editing the documentation
=========================

After modifying some files in the Sage tutorial
(:sage_root:`src/doc/en/tutorial/`), you will want to visualize the result. In
order to build an **HTML** version of this document, type

.. code-block:: console

    $ meson compile -C builddir doc-html-other-en-tutorial

Here and below, replace ``builddir`` with the Meson build directory.  An
editable pip installation uses ``build/cpXY``, where ``X`` and ``Y`` are the
major and minor Python versions; for example, Python 3.12 uses ``build/cp312``.
A build using Sage's top-level Makefile uses ``build/sage-distro``.  See
:ref:`section-meson-build-directory` for details.  You can now open
:file:`builddir/src/doc/html/en/tutorial/index.html` in your web browser.

- Do you want to **add a new file** to the documentation? :ref:`Click here
  <section-add-file>`.

- For more detailed information on the documentation build targets, see
  :ref:`section-building-manuals`.

**Run doctests:** All files must pass tests. After modifying a document
(e.g. ``tutorial``), you can run tests with the following command (see
:ref:`chapter-testing`)

.. code-block:: console

    $ sage -tp SAGE_ROOT/src/doc/en/tutorial/

**Reference manual:** as this manual is mostly generated from Sage's source
code, you will need to build Sage in order to see the changes you made to some
function's documentation.  Type

.. code-block:: console

    $ meson compile -C builddir
    $ meson compile -C builddir doc-html-reference-reference_top

.. _chapter-sage_manuals_links:

Hyperlinks
==========

The documentation can contain links toward modules, classes, or methods, e.g.:

.. CODE-BLOCK:: rest

    :mod:`link to a module <sage.module_name>`
    :mod:`sage.module_name` (here the link's text is the module's name)

For links toward classes, methods, or functions, replace ``:mod:`` by
``:class:``, ``:meth:``, or ``:func:``, respectively.  See Sphinx' documentation
on `cross-referencing Python objects
<https://www.sphinx-doc.org/en/master/usage/domains/python.html#cross-referencing-python-objects>`_
and for the general syntax of
`roles <https://www.sphinx-doc.org/en/master/usage/restructuredtext/roles.html>`_.

**Short links:** the link ``:func:`~sage.mod1.mod2.mod3.func1``` is equivalent
to ``:func:`func1 <sage.mod1.mod2.mod3.func1>```: the function's name will be
used as the link name, instead of its full path.

**Local names:** links between methods of the same class do not need to be
absolute. If you are documenting ``method_one``, you can write
``:meth:`method_two```.

**Intersphinx references:** in the same way, you can refer to the modules, classes,
methods, functions of the Python standard library and of several Python packages
used by SageMath; see the `Intersphinx documentation
<https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html>`_
for details. Likewise, you can refer to the C functions of the
:ref:`FLINT <spkg_flint>` library; see `Sphinx' documentation on
cross-referencing C constructs
<https://www.sphinx-doc.org/en/master/usage/domains/c.html#cross-referencing-c-constructs>`_
for more information.

.. LIST-TABLE::
   :widths: 4 7 5
   :header-rows: 0

   * - Python
     - ``:exc:`ValueError```
     - :exc:`ValueError`
   * - :ref:`CVXOPT <spkg_cvxopt>`
     - ``:func:`cvxopt.solvers.socp```
     - :func:`cvxopt.solvers.socp`
   * - :ref:`CVXpy <spkg_cvxpy>`
     - ``:class:`~cvxpy.atoms.log_det.log_det```
     - :class:`~cvxpy.atoms.log_det.log_det`
   * - :ref:`cypari2 <spkg_cypari>`
     - ``:class:`cypari2.gen.Gen```
     - :class:`cypari2.gen.Gen`
   * - :ref:`cysignals <spkg_cysignals>`
     - ``:envvar:`CYSIGNALS_CRASH_DAYS```
     - :envvar:`CYSIGNALS_CRASH_DAYS`
   * - :ref:`FLINT <spkg_flint>`
     - ``:c:func:`arith_bell_number```
     - :c:func:`arith_bell_number`
   * - :ref:`gmpy2 <spkg_gmpy2>`
     - ``:func:`gmpy2.gamma_inc```
     - :func:`gmpy2.gamma_inc`
   * - :ref:`ipywidgets <spkg_ipywidgets>`
     - ``:mod:`~ipywidgets.widgets.widget_date```
     - :mod:`~ipywidgets.widgets.widget_date`
   * - :ref:`Matplotlib <spkg_matplotlib>`
     - ``:mod:`matplotlib.bezier```
     - :mod:`matplotlib.bezier`
   * - :ref:`mpmath <spkg_mpmath>`
     - ``:attr:`mpmath.mp.khinchin```
     - :attr:`mpmath.mp.khinchin`
   * - :ref:`NetworkX <spkg_networkx>`
     - ``:attr:`~networkx.DiGraph.out_degree```
     - :attr:`~networkx.DiGraph.out_degree`
   * - :ref:`NumPy <spkg_numpy>`
     - ``:data:`numpy.NAN```
     - :data:`numpy.NAN`
   * - :ref:`pplpy <spkg_pplpy>`
     - ``:mod:`ppl.polyhedron```
     - :mod:`ppl.polyhedron`
   * - :ref:`rpy2 <spkg_rpy2>`
     - ``:class:`~rpy2.robjects.vectors.DataFrame```
     - :class:`~rpy2.robjects.vectors.DataFrame`
   * - :ref:`SciPy <spkg_scipy>`
     - ``:data:`scipy.special.huber```
     - :data:`scipy.special.huber`
   * - :ref:`SymPy <spkg_sympy>`
     - ``:class:`~sympy.diffgeom.WedgeProduct```
     - :class:`~sympy.diffgeom.WedgeProduct`

To see the available cross references in any of these libraries, you can use the command
``./sage -python -m sphinx.ext.intersphinx src/doc/common/_vendor/numpy.inv``.

**Global namespace:** if an object (e.g. ``integral``) is automatically imported
by Sage, you can link toward it without specifying its full path:

.. CODE-BLOCK:: rest

    :func:`A link toward the integral function <integral>`

**Sage-specific roles:** Sage defines several specific *roles*:

.. LIST-TABLE::
   :widths: 4 4 4
   :header-rows: 0

   * - GitHub issue
     - ``:issue:`17596```
     - :issue:`17596`

   * - Sage repository file or directory
     - ``:sage_root:`src/doc/en```
     - :sage_root:`src/doc/en`

   * - Wikipedia
     - ``:wikipedia:`Sage_(mathematics_software)```
     - :wikipedia:`Sage_(mathematics_software)`

   * - arXiv
     - ``:arxiv:`1202.1506```
     - :arxiv:`1202.1506`

   * - On-Line Encyclopedia of Integer Sequences
     - ``:oeis:`A000081```
     - :oeis:`A000081`

   * - Digital Object Identifier
     - ``:doi:`10.2752/175303708X390473```
     - :doi:`10.2752/175303708X390473`

   * - MathSciNet
     - ``:mathscinet:`MR0100971```
     - :mathscinet:`MR0100971`

   * - :ref:`ECL <spkg_ecl>`
     - ``:ecl:`Manipulating-Lisp-objects```
     - :ecl:`Manipulating-Lisp-objects`

   * -
     - ``:common_lisp:`RENAME-PACKAGE <f_rn_pkg>```
     - :common_lisp:`RENAME-PACKAGE <f_rn_pkg>`

   * - :ref:`GAP <spkg_gap>`
     - ``:gap:`Groups <chap39>```
     - :gap:`Groups <chap39>`

   * -
     - ``:gap_package:`GAP package QuaGroup <quagroup/doc/chap0_mj.html>```
     - :gap_package:`GAP package QuaGroup <quagroup/doc/chap0_mj.html>`

   * - :ref:`Giac <spkg_giac>`
     - ``:giac_cascmd:`gbasis <node280>```
     - :giac_cascmd:`gbasis <node280>`

   * -
     - ``:giac_us:`Unary-functions```
     - :giac_us:`Unary-functions`

   * - :ref:`Maxima <spkg_maxima>`
     - ``:maxima:`struve_h <index-struve_005fh>```
     - :maxima:`struve_h <index-struve_005fh>`

   * - :ref:`Meson <spkg_meson>`
     - ``:meson:`install_subdir <Reference-manual_functions.html#install_subdir>```
     - :meson:`install_subdir <Reference-manual_functions.html#install_subdir>`

   * - :ref:`Pari <spkg_pari>`
     - ``:pari:`lfungenus2```
     - :pari:`lfungenus2`

   * - :ref:`polymake <spkg_polymake>`
     - ``:polymake:`matroid```
     - :polymake:`matroid`

   * - :ref:`PPL <spkg_ppl>`
     - ``:ppl:`Linear_Expression <classParma__Polyhedra__Library_1_1 Linear__Expression>```
     - :ppl:`Linear_Expression <classParma__Polyhedra__Library_1_1Linear__Expression>`

   * - :ref:`QEPCAD <spkg_qepcad>`
     - ``:qepcad:`QEPCAD: Entering formulas <user/EnterForm>```
     - :qepcad:`QEPCAD: Entering formulas <user/EnterForm>`

   * - :ref:`SCIP <spkg_scip>`
     - ``:scip:`SCIPsolve <group__PublicSolveMethods>```
     - :scip:`SCIPsolve <group__PublicSolveMethods>`

   * - :ref:`Singular <spkg_singular>`
     - ``:singular:`stdfglm <sing_358>```
     - :singular:`stdfglm <sing_358>`

   * - :ref:`SoPlex <spkg_soplex>`
     - ``:soplex:`soplex::LinSolverRational <classsoplex_1_1SLinSolverRational>```
     - :soplex:`soplex::LinSolverRational <classsoplex_1_1SLinSolverRational>`

**http links:** copy/pasting a http link in the documentation works. If you want
a specific link name, use ```link name <http://www.example.com>`_``

**Anonymous hyperlinks:** Using a single underscore creates an *explicit target
name* ``"link name"`` which needs to be unique in the current page. Using the
same target name twice in the same page creates an error while building the
documentation saying ``WARNING: Duplicate explicit target name: ...``. To
avoid this issue, one can change the target names to be all different or
another option is to use `anonymous hyperlinks
<https://stackoverflow.com/questions/27420317/>`__ with two underscores, as in
``see `this page <http://www.example.com>`__ or `this page
<http://www.example2.com>`__``.

**Broken links:** Sphinx can report broken links. See
:ref:`section-building-manuals`.

.. _section-add-file:

Adding a new file
=================

If you added a new file to Sage (e.g. ``sage/matroids/my_algorithm.py``) and you
want its content to appear in the reference manual, you have to add its name to
the file :sage_root:`src/doc/en/reference/matroids/index.rst`. Replace
'matroids' with whatever fits your case.

**The combinat/ folder:** if your new file belongs to a subdirectory of combinat/ the
procedure is different:

* Add your file to the index stored in the ``__init__.py`` file located in the
  directory that contains your file.

* Add your file to the index contained in
  :sage_root:`src/doc/en/reference/combinat/module_list.rst`.

.. _section-building-manuals:

Building the manuals
====================

*(Do you want to edit the documentation?* :ref:`Click here
<section-manuals-edit>`)

All Sage manuals are built through Meson targets.  These targets invoke the
documentation builder in :sage_root:`src/sage_docbuild`, which is a wrapper
around Sphinx.  The pip-installed ``sage`` command does not provide the legacy
``--docbuild`` option, so use the Meson targets instead.

In the commands below, replace ``builddir`` with the Meson build directory.
For an editable pip installation this is ``build/cpXY``, such as
``build/cp312`` for Python 3.12.  A build using Sage's top-level Makefile uses
``build/sage-distro``; see :ref:`section-meson-build-directory`.  To build all
HTML or PDF manuals, use

.. code-block:: console

    $ meson compile -C builddir doc-html
    $ meson compile -C builddir doc-pdf

Building PDF manuals requires a working LaTeX installation.  The generated
manuals are written below :file:`builddir/src/doc/html` and
:file:`builddir/src/doc/pdf`, respectively.

Meson also provides targets for individual manuals.  To list all configured
targets, use

.. code-block:: console

    $ meson introspect builddir --targets

**Broken links:** in order to build the documentation while reporting the broken
links that it contains, pass ``--warn-links`` through
``SAGE_DOCBUILD_OPTS``. Note that Sphinx will not rebuild a document that has
not been updated, and thus not report its broken links.

.. code-block:: console

    $ SAGE_DOCBUILD_OPTS="--warn-links" meson compile -C builddir doc-html

.. _section-manuals-names:

Document names
--------------

Documentation sources outside the reference manual have the form

.. CODE-BLOCK:: text

    lang/name

where ``lang`` is a two-letter language code, and ``name`` is the
descriptive name of the document.  The corresponding HTML target is named
``doc-html-other-lang-name``.  For example, to build the English or French
tutorial, use

.. code-block:: console

    $ meson compile -C builddir doc-html-other-en-tutorial
    $ meson compile -C builddir doc-html-other-fr-tutorial

Reference manual targets omit the language and ``other`` components.  For
example, to build the combinatorics part of the reference manual, use

.. code-block:: console

    $ meson compile -C builddir doc-html-reference-combinat

Replace ``html`` with ``pdf`` in these target names to build PDF output.


Syntax highlighting Cython code
===============================

If you want to write :ref:`Cython <chapter-cython>` code in a ReST file, precede
the code block by ``.. CODE-BLOCK:: cython`` instead of the usual ``::``. Enable
syntax-highlighting in a whole file with ``.. HIGHLIGHT:: cython``. Example:

.. CODE-BLOCK:: cython

    cdef extern from "descrobject.h":
        ctypedef struct PyMethodDef:
            void *ml_meth
        ctypedef struct PyMethodDescrObject:
            PyMethodDef *d_method
        void* PyCFunction_GET_FUNCTION(object)
        bint PyCFunction_Check(object)
