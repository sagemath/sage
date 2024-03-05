=========================================================================================
 Sage: Open Source Mathematics Software: Script for maintaining a SageMath-based project
=========================================================================================

About SageMath
--------------

   "Creating a Viable Open Source Alternative to
    Magma, Maple, Mathematica, and MATLAB"

   Copyright (C) 2005-2024 The Sage Development Team

   https://www.sagemath.org

SageMath fully supports all major Linux distributions, recent versions of macOS, and Windows (using Windows Subsystem for Linux).

The traditional and recommended way to install SageMath is from source via Sage-the-distribution (https://www.sagemath.org/download-source.html).  Sage-the-distribution first builds a large number of open source packages from source (unless it finds suitable versions installed in the system) and then installs the Sage Library (sagelib, implemented in Python and Cython).


About this pip-installable source distribution
----------------------------------------------

Creating a user project
~~~~~~~~~~~~~~~~~~~~~~~

::

   $ sage-project create PROJECT-DIRECTORY

This creates configuration files:

- ``environment*.yml`` for local use with conda-forge in Linux, macOS
- ``.devcontainer/downstream-*/`` for use in dev container on Linux, macOS, Windows

It can also be invoked as follows::

   $ pipx run cookiecutter gh:sagemath/sage --directory="pkgs/sage-project/sage_project/user-project-template"

See https://cookiecutter.readthedocs.io/en/latest/README.html for available options.


Creating a pip-installable downstream package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  $ sage-project create --downstream-package PROJECT-DIRECTORY

Additionally creates:

- ``.github/workflows/``


Adding Sage CI portability/integration testing infrastructure to an upstream project
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   $ sage-project create --upstream-package PROJECT-DIRECTORY

Creates:

- ``.github/workflows/ci-sage.yml``
- ``.devcontainer/portability-*``
- ``.devcontainer/tox-docker-in-docker``
