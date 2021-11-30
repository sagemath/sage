.. highlight:: bash

.. _chapter-sage_manuals:

================
The Sage Manuals
================

Sage's manuals are written in `ReST <http://docutils.sourceforge.net/rst.html>`_
(reStructuredText), and generated with the software `Sphinx
<http://sphinx.pocoo.org>`_:

.. LIST-TABLE::
   :widths: 4 12
   :header-rows: 1

   * - Name
     - Files

   * - `Tutorial <../tutorial/index.html>`_
     - ``SAGE_ROOT/src/doc/en/tutorial``

   * - `Developer's guide <../developer/index.html>`_
     - ``SAGE_ROOT/src/doc/en/developer``

   * - `Constructions <../constructions/index.html>`_
     - ``SAGE_ROOT/src/doc/en/constructions``

   * - `Installation guide <../installation/index.html>`_
     - ``SAGE_ROOT/src/doc/en/installation``

   * - `Reference manual <../reference/index.html>`_
     - ``SAGE_ROOT/src/doc/en/reference``
       (most of it is generated from the
       source code)

- Additionally, more specialized manuals can be found under
  ``SAGE_ROOT/src/doc/en``.

- Some documents have been **translated** into other languages. In order to
  access them, change en/ into fr/,es/, de/... See :ref:`section-manuals-names`.

.. _section-manuals-edit:

Editing the documentation
=========================

After modifying some files in the Sage tutorial
(``SAGE_ROOT/src/doc/en/tutorial/``), you will want to visualize the result. In
order to build a **html** version of this document, type::

    sage --docbuild tutorial html

You can now open ``SAGE_ROOT/local/share/doc/sage/html/en/tutorial/index.html`` in
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

For links toward classes, methods, or function, replace **:mod:** by
**:class:**, **:meth:** or **func:** respectively.  See `Sphinx' documentation
<http://sphinx.pocoo.org/markup/inline.html>`_.

**Short links:** the link ``:func:`~sage.mod1.mod2.mod3.func1``` is equivalent
to ``:func:`func1 <sage.mod1.mod2.mod3.func1>```: the function's name will be
used as the link name, instead of its full path.

**Local names:** links between methods of the same class do not need to be
absolute. If you are documenting ``method_one``, you can write
``:meth:`method_two```.

**Global namespace:** if an object (e.g. ``integral``) is automatically imported
by Sage, you can link toward it without specifying its full path:

.. CODE-BLOCK:: rest

    :func:`A link toward the integral function <integral>`

**Sage-specific roles:** Sage defines several specific *roles*:

.. LIST-TABLE::
   :widths: 4 4 4
   :header-rows: 0

   * - Trac server
     - ``:trac:`17596```
     - :trac:`17596`

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

Adding a New File
=================

If you added a new file to Sage (e.g. ``sage/matroids/my_algorithm.py``) and you
want its content to appear in the reference manual, you have to add its name to
the file ``SAGE_ROOT/src/doc/en/reference/matroids/index.rst``. Replace
'matroids' with whatever fits your case.

**The combinat/ folder:** if your new file belongs to a subdirectory of combinat/ the
procedure is different:

* Add your file to the index stored in the ``__init__.py`` file located in the
  directory that contains your file.

* Add your file to the index contained in
  ``SAGE_ROOT/src/doc/en/reference/combinat/module_list.rst``.

.. _section-building-manuals:

Building the Manuals
====================

*(Do you want to edit the documentation?* :ref:`Click here
<section-manuals-edit>`)

All of the Sage manuals are built using the ``sage --docbuild``
script.  The content of the ``sage --docbuild`` script is defined in
``SAGE_ROOT/src/sage_docbuild/__init__.py``.  It is a thin wrapper around
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
used in Sage. See `<http://sphinx.pocoo.org/builders.html>`_.

**Broken links:** in order to build the documentation while reporting the broken
links that it contains, use the ``--warn-links`` flag. Note that Sphinx will not
rebuild a document that has not been updated, and thus not report its broken
links::

        sage --docbuild --warn-links reference html

.. _section-manuals-names:

Document Names
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


Syntax Highlighting Cython Code
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
