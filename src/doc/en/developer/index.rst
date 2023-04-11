.. _developers-guide:

======================================
Welcome to the Sage Developer's Guide!
======================================

.. WARNING::

    **Sage development moved to GitHub in February 2023.** After the transition,
    some parts of this guide (especially those related with `the Sage Trac
    server <https://trac.sagemath.org>`_) became obsolete and need to be
    updated according to the new workflow on GitHub. See our `transition guide
    from Trac to GitHub
    <https://github.com/sagemath/trac-to-github/blob/master/docs/Migration-Trac-to-Github.md>`_
    for the preliminary version of the workflow.

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

- **Trac server:** all changes go through `the Sage Trac server
  <https://trac.sagemath.org>`_ at some point. It contains bug reports, upgrade
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

  As an easy way to get started, you can run and edit Sage's code and contribute
  your changes using `Gitpod <https://www.gitpod.io>`_,
  a free online development environment based on VS Code.
  It will launch a pre-made workspace with all dependencies and tools installed
  so that you can start contributing straight away.
  Start by `going to Gitpod <https://gitpod.io/#https://github.com/sagemath/sage>`_,
  and read :ref:`our Gitpod guidelines <section-gitpod>` to learn more.

- **Conventions:** read our :ref:`conventions and guidelines
  <section-writing-code-for-sage>` for code and documentation.

  For everything related to manuals, tutorials, and languages, :ref:`click here
  <chapter-sage_manuals>`.

- **Git (revision control):** To share changes with the Sage community, you will
  need to learn about revision control; we use the software Git for this
  purpose.

  - :ref:`How to install it? <section-git-install>`
  - :ref:`How to configure it for use with Trac? <section-git-setup-name>`
  - :ref:`Here is <chapter-walkthrough>` an overview of our development flow.

Git and Trac for Sage development
=================================

First Steps with Git
--------------------

Sage uses git for version control.

.. toctree::
   :maxdepth: 3

   git_setup
   walk_through

.. _section-git-tricks-and-tips:

Using Git with Trac
-------------------

To contribute back your changes to Sage source code to the project,
you will need a ticket on the
`Sage trac server <http://trac.sagemath.org>`_.

.. toctree::
   :maxdepth: 2

   trac
   manual_git
   git_background
   advanced_git
   workflows
   git_trac


.. _section-writing-code-for-sage:

Writing Code for Sage
=====================

.. toctree::
   :maxdepth: 3

   workspace
   coding_basics
   reviewer_checklist

Running Sage's tests
--------------------

.. toctree::
   :maxdepth: 3

   doctesting

Testing on multiple platforms
-----------------------------

.. toctree::
   :maxdepth: 3

   portability_testing

Additional development and testing tools
----------------------------------------

.. toctree::
   :maxdepth: 3

   tools

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

Packaging the Sage Library
--------------------------

.. toctree::
   :maxdepth: 3

   packaging_sage_library

Packaging Third-Party Code
--------------------------

.. toctree::
   :maxdepth: 3

   packaging


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License <http://creativecommons.org/licenses/by-sa/3.0/>`_.
