======================================
Welcome to the Sage Developer's Guide!
======================================

Everybody who uses Sage is encouraged to contribute something back to Sage at
some point. You could:

* Add examples to the documentation
* Find bugs or typos
* Fix a bug
* Implement a new function
* Contribute a useful tutorial for a mathematical topic
* Translate an existing document to a new language
* Create a new class, create a fast new C library, etc.

This document tells you what you need to know to do all the above, from
reporting bugs to modifying and extending Sage and its documentation.  We also
discuss how to share your new and modified code with other Sage users around the
globe.

Here are brief overviews of each part; for more details, see the extended table
of contents below.  No matter where you start, good luck and welcome to Sage
development!

- **Trac server:** all changes go through the `the Sage Trac server
  <http://trac.sagemath.org>`_ at some point. It contains bug reports, upgrade
  requests, changes in progress, and those already part of Sage
  today. :ref:`Click here <chapter-sage-trac>` for more information.

  Importantly, you will need to :ref:`create a trac account
  <section-trac-account>` in order to contribute.

- **Source code:** You need your own copy of Sage's source code to change it.
  `Go there <http://doc.sagemath.org/html/en/installation/source.html>`_ to get it
  and for instructions to build it.

  If you have never worked on software before, pay close attention to the
  `prerequisites to compile
  <http://doc.sagemath.org/html/en/installation/source.html#prerequisites>`_ on your
  system.

- **Conventions:** read our :ref:`conventions and guidelines
  <section-writing-code-for-sage>` for code and documentation.

  For everything related to manuals, tutorials, and languages, :ref:`click here
  <chapter-sage_manuals>`.

- **Git (revision control):** To share changes with the Sage community, you will
  need to learn about revision control; we use the software Git for this
  purpose.

  - :ref:`Here is <chapter-walkthrough>` an overview of our development flow.
  - :ref:`Unfamiliar with Git or revision control? <chapter-git_trac>`
  - :ref:`How to install it? <section-git-install>`
  - :ref:`How to configure it for use with Trac? <section-git-setup-name>`

Git for Sage development
========================

First Steps with Git
--------------------

Sage uses git for version control.

.. toctree::
   :maxdepth: 3

   git_setup
   walk_through

The git-trac command
--------------------

Putting your local changes on a Trac ticket.

.. toctree::
   :maxdepth: 2

   git_trac

.. _section-git-tricks-and-tips:

Git Tricks & Tips
-----------------

When ``git trac`` is not enough.

.. toctree::
   :maxdepth: 2

   manual_git
   git_background
   advanced_git
   workflows

Sage Trac and tickets
=====================

All changes to Sage source code require a ticket on the
`Sage trac server <http://trac.sagemath.org>`_.

.. toctree::
   :maxdepth: 2

   trac


.. _section-writing-code-for-sage:

Writing Code for Sage
=====================

.. toctree::
   :maxdepth: 3

   coding_basics
   reviewer_checklist

Running Sage's tests
--------------------

.. toctree::
   :maxdepth: 3

   doctesting

Contributing to Manuals and Tutorials
-------------------------------------

.. toctree::
   :maxdepth: 3

   sage_manuals

Sage Coding Details
-------------------

.. toctree::
   :maxdepth: 3

   coding_in_python
   coding_in_cython
   coding_in_other

Packaging Third-Party Code
--------------------------

.. toctree::
   :maxdepth: 3

   packaging
   packaging_old_spkgs

Sage Notebook Developer Guide
=============================

.. toctree::
   :maxdepth: 3

   sagenb/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License <http://creativecommons.org/licenses/by-sa/3.0/>`_.
