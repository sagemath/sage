.. _developers-guide:

===============================
Welcome to Sage Developer Guide
===============================

Everybody who uses Sage is encouraged to contribute something back to Sage at
some point. You could:

* Add examples to the documentation
* Find bugs or typos
* Fix a bug
* Implement a new function or create a new class
* Contribute a useful tutorial for a mathematical topic
* Translate an existing document to a new language
* Upgrade a package, create a fast new C library, etc.

This document tells you what you need to know to do all the above. We also
discuss how to share your new and modified code with other Sage users around
the globe.

To begin with, you need of course your own copy of Sage source code to change
it. Use our `Installation guide
<http://doc.sagemath.org/html/en/installation/source.html>`_ to get the source
code and build Sage from source. If you have never worked on software before,
pay close attention to the prerequisites to build on your platform.

Now here is a brief overview of this guide.

- :ref:`section-first-steps`: To share changes with the Sage community, you
  need to learn about revision control. We use the software Git for this
  purpose. Here we walk you through from setting up Git on your platform
  and to preparing a local branch to share with all Sage users.

  .. NOTE::

    As an easy way to get started, you can run and edit Sage's code and contribute
    your changes using `Gitpod <https://www.gitpod.io>`_, a free online development
    environment based on VS Code.  It will launch a pre-made workspace with all
    dependencies and tools installed so that you can start contributing straight
    away.  Start by `going to Gitpod
    <https://gitpod.io/#https://github.com/sagemath/sage>`_, and read :ref:`our
    Gitpod guidelines <section-gitpod>` to learn more.

- :ref:`section-development-on-github`: All changes go through `the Sage
  repository on GitHub <https://github.com/sagemath/sage>`_ at some point. It contains
  bug reports, enhancement proposals, changes in progress, and indeed all the
  history of Sage today. You have to be familiar with it to be involved in Sage
  development.

- :ref:`section-git-tricks-and-tips`: Here we give an in-depth guide for
  working with Git for Sage development. Read this when you need help on Git in
  a tricky situation such as merge conflict.

- :ref:`section-writing-code-for-sage`: This is a guide on conventions in
  writing code and documentation. A beginning developer should read this to be
  a good developer. As conventions evolve over time, also experienced Sage
  contributors may want to review this chapter once in a while.

- :ref:`section-testing-sage`: We value testing Sage highest. Every change of
  Sage source code has a risk to break Sage, and must be tested before being
  merged.  This part explains our various tools to help test Sage.

- :ref:`section-updating-documentation`: All features of Sage are documented in
  our manuals. This part explains the technical aspect of updating Sage
  documentation.

- :ref:`section-more-on-coding`: When you need to know the technical details of
  Sage for deep coding, read this.

- :ref:`section-packaging`: Sage is composed of many third-party packages and
  its own distribution packages. This part is for advanced developers.

For more details, see the table of contents below. No matter where you start,
good luck and welcome to Sage development!


Table of Contents
=================

.. _section-first-steps:

First Steps
-----------

.. toctree::
   :maxdepth: 2

   walkthrough
   git_setup


.. _section-development-on-github:

Working on GitHub
-----------------

.. toctree::
   :maxdepth: 2

   github
   workflows
   review


.. _section-git-tricks-and-tips:

Working with Git
----------------

.. toctree::
   :maxdepth: 2

   git_basic
   git_advanced
   git_background


.. _section-writing-code-for-sage:

Writing Code for Sage
---------------------

.. toctree::
   :maxdepth: 2

   workspace
   coding_basics


.. _section-testing-sage:

Testing Sage
------------

.. toctree::
   :maxdepth: 2

   doctesting
   portability_testing
   tools


.. _section-updating-documentation:

Updating Sage Documentation
---------------------------

.. toctree::
   :maxdepth: 2

   sage_manuals


.. _section-more-on-coding:

More on Coding for Sage
-----------------------

.. toctree::
   :maxdepth: 2

   coding_in_python
   coding_in_cython
   coding_in_other


.. _section-packaging:

Packaging
---------

.. toctree::
   :maxdepth: 2

   packaging
   downstream
   packaging_sage_library


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License <http://creativecommons.org/licenses/by-sa/3.0/>`_.
