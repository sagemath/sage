.. _section-command-line:

The Sage Command Line
=====================

The Sage Read-Eval-Print-Loop (REPL) is based on IPython. In this
document, you'll find how the IPython integration works. You should
also be familiar with the documentation for IPython.

For more details about using the Sage command line, see :ref:`the Sage
tutorial <../../tutorial/index.html>`.

Running Sage
------------

.. toctree::
   :maxdepth: 1

   options
   startup
   environ
   sage/misc/trace

Preparsing
----------

Sage commands are "preparsed" to valid Python syntax. This allows
for example to support the ``R.<x> = QQ[]`` syntax.

.. toctree::
   :maxdepth: 1

   sage/repl/preparse


Loading and attaching files
---------------------------

Sage or Python files can be loaded (similar to Python's ``execfile``)
in a Sage session. Attaching is similar, except that the attached file
is reloaded whenever it is changed.

.. toctree::
   :maxdepth: 1

   sage/repl/load
   sage/repl/attach


Pretty Printing
---------------

In addition to making input nicer, we also modify how results are
printed. This again builds on how IPython formats output. Technically,
this works using a modified displayhook in Python.

.. toctree::
   :maxdepth: 1

   sage/repl/display/formatter
   sage/repl/display/pretty_print
   sage/repl/display/fancy_repr
   sage/repl/display/util

Display Backend Infrastructure
------------------------------

.. toctree::
   :maxdepth: 1

   sage/repl/rich_output/display_manager
   sage/repl/rich_output/preferences
   sage/repl/rich_output/pretty_print
   sage/repl/rich_output/buffer
   sage/repl/rich_output/output_basic
   sage/repl/rich_output/output_graphics
   sage/repl/rich_output/output_graphics3d
   sage/repl/rich_output/output_video
   sage/repl/rich_output/output_catalog

   sage/repl/rich_output/backend_base
   sage/repl/rich_output/test_backend
   sage/repl/rich_output/backend_doctest
   sage/repl/rich_output/backend_ipython

Miscellaneous
-------------

.. toctree::
   :maxdepth: 1

   sage/repl/interpreter
   sage/repl/ipython_extension
   sage/repl/interface_magic
   sage/repl/ipython_kernel/interact
   sage/repl/ipython_kernel/widgets
   sage/repl/ipython_kernel/install
   sage/repl/ipython_kernel/kernel
   sage/repl/ipython_tests

   sage/repl/display/jsmol_iframe
   sage/repl/image
   sage/repl/inputhook

.. include:: ../footer.txt
