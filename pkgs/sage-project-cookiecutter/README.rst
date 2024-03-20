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

   $ sage-project-cookiecutter create PROJECT-DIRECTORY

This creates configuration files:

- ``environment*.yml`` for local use with conda-forge in Linux, macOS
- ``.devcontainer/downstream-*/`` for use in dev container on Linux, macOS, Windows

It can also be invoked as follows::

   $ pipx run cruft create https://github.com/mkoeppe/sage --checkout sagemath-environment-cookiecutter \
       --directory="pkgs/sage-project-cookiecutter/sage_project_cookiecutter/user-project-template"

TODO: Update URL, remove --checkout ... before merging.

See https://cruft.github.io/cruft/ for available options.


Creating a pip-installable downstream package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  $ sage-project-cookiecutter create --downstream-package PROJECT-DIRECTORY

Additionally creates:

- ``.github/workflows/``

It can also be invoked as follows::

   $ pipx run cruft create https://github.com/mkoeppe/sage --checkout sagemath-environment-cookiecutter \
       --directory="pkgs/sage-project-cookiecutter/sage_project_cookiecutter/downstream-package-template"


Adding Sage CI portability/integration testing infrastructure to an upstream project
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   $ sage-project-cookiecutter create --upstream-package PROJECT-DIRECTORY

Creates in the existing ``PROJECT-DIRECTORY``:

- ``.github/workflows/ci-sage.yml``
- ``.devcontainer/portability-*``
- ``.devcontainer/tox-docker-in-docker``

It can also be invoked as follows::

   [alice@localhost PROJECT-DIRECTORY]$ (cd .. && pipx run cruft create \
       https://github.com/mkoeppe/sage \
       --checkout sagemath-environment-cookiecutter \
       --directory="pkgs/sage-project-cookiecutter/sage_project_cookiecutter/upstream-package-template" \
       --overwrite-if-exists)
   [1/1] Name of the project (directory name to create) (my-sage-project): PROJECT-DIRECTORY


Creating a pip-installable upstream package of the SageMath organization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   $ sage-project-cookiecutter create --sagemath-upstream-package PROJECT-DIRECTORY

Additionally creates:

- ``CODE_OF_CONDUCT.md``
- ``CONTRIBUTING.md``

It can also be invoked as follows::

   $ pipx run cruft create https://github.com/mkoeppe/sage --checkout sagemath-environment-cookiecutter \
       --directory="pkgs/sage-project-cookiecutter/sage_project_cookiecutter/sagemath-upstream-package-template"


Updating a project
~~~~~~~~~~~~~~~~~~

::

   [alice@localhost PROJECT-DIRECTORY]$ pipx run cruft update \
       --checkout sagemath-environment-cookiecutter
