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
order to build a **html** version of this document, type::

    sage --docbuild tutorial html

You can now open :file:`SAGE_ROOT/local/share/doc/sage/html/en/tutorial/index.html` in
your web browser.

- Do you want to **add a new file** to the documentation? :ref:`Click here
  <section-add-file>`.

- For more detailed information on the ``--docbuild`` command, see
  :ref:`section-building-manuals`.

**Run doctests:** All files must pass tests. After modifying a document
(e.g. ``tutorial``), you can run tests with the following command (see
:ref:`chapter-testing`)::

    sage -tp SAGE_ROOT/src/doc/en/tutorial/

**Reference manual:** as this manual is mostly generated from Sage's source
code, you will need to build Sage in order to see the changes you made to some
function's documentation.  Type::

    sage -b && sage --docbuild reference html

.. _chapter-sage_manuals_links:

Hyperlinks
==========

The documentation can contain links toward modules, classes, or methods, e.g.::

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

   * - :ref:`CMR <spkg_cmr>`
     - ``:cmr:`GraphNode <structGraphNode>```
     - :cmr:`GraphNode <structGraphNode>`

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

.. _section-documentation-conditional:

Making portions of the reference manual conditional on optional features
========================================================================

For every dynamically detectable feature such as :class:`graphviz
<~sage.features.graphviz.Graphviz>` or :class:`sage.symbolic
<sage.features.sagemath.sage__symbolic>` (see :mod:`sage.features`),
Sage defines a Sphinx tag that can be used with the `Sphinx
directive ".. ONLY::"
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#tags>`_.
Because Sphinx tags have to use Python identifier syntax, Sage uses
the format ``feature_``, followed by the feature name where dots are
replaced by underscores. Hence, conditionalizing on the features of
the previous examples would look as follows:

.. CODE-BLOCK:: rest

  .. ONLY:: feature_graphviz

and:

.. CODE-BLOCK:: rest

  .. ONLY:: feature_sage_symbolic

.. _section-building-manuals:

Building the manuals
====================

*(Do you want to edit the documentation?* :ref:`Click here
<section-manuals-edit>`)

All of the Sage manuals are built using the ``sage --docbuild``
script.  The content of the ``sage --docbuild`` script is defined in
:sage_root:`src/sage_docbuild/__init__.py`.  It is a thin wrapper around
the ``sphinx-build`` script which does all of the real work.  It is
designed to be a replacement for the default Makefiles generated by
the ``sphinx-quickstart`` script.  The general form of the command
is::

    sage --docbuild <document-name> <format>

For example::

    sage --docbuild reference html

Two **help** commands which give plenty of documentation for the ``sage
--docbuild`` script::

    sage --docbuild -h # short help message
    sage --docbuild -H # a more comprehensive one

**Output formats:** All output formats supported by Sphinx (e.g. pdf) can be
used in Sage. See `<http://www.sphinx-doc.org/builders.html>`_.

**Broken links:** in order to build the documentation while reporting the broken
links that it contains, use the ``--warn-links`` flag. Note that Sphinx will not
rebuild a document that has not been updated, and thus not report its broken
links::

        sage --docbuild --warn-links reference html

.. _section-manuals-names:

Document names
--------------

The ``<document-name>`` has the form:

.. CODE-BLOCK:: text

    lang/name

where ``lang`` is a two-letter language code, and ``name`` is the
descriptive name of the document.  If the language is not specified,
then it defaults to English (``en``).  The following two commands do
the exact same thing::

    sage --docbuild tutorial html
    sage --docbuild en/tutorial html

To specify the French version of the tutorial, you would simply run::

    sage --docbuild fr/tutorial html


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
