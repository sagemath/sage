.. _spkg:

Packages and Features
=====================

Standard Packages
-----------------

The Sage distribution includes most programs and libraries on which
Sage depends.  It installs them automatically if it does not find
equivalent system packages.

.. envvar:: FC

   Fortran compiler command passed to ``./configure``.

.. include:: index_standard.rst


Optional Packages
-----------------

For additional functionality, you can install some of the following
optional packages.

.. include:: index_optional.rst


Installing Optional Packages
----------------------------

To install an optional package, you can use Sage's package management system.

Basic Installation
~~~~~~~~~~~~~~~~~~

To install an optional package, use the following command in a system
terminal (not from the Sage prompt):

.. code-block:: console

    sage -i <package_name>

For example, to install the optional package `bliss`:

.. code-block:: console

    sage -i bliss

.. note::

    The ``sage -i`` command is part of Sage-the-distribution. It may not
    work in Conda or other packaged installations. In that case, use the
    package manager of your installation instead.

Using the Package List
~~~~~~~~~~~~~~~~~~~~~~

To see a list of all available optional packages, run in a system terminal:

.. code-block:: console

    sage --optional

To list experimental packages instead:

.. code-block:: console

    sage --experimental

For more detailed information on listing packages, including how to list
pip-installed packages, see the documentation for the
:mod:`sage.misc.package` module.

Alternatively, inside Sage you can use:

.. code-block:: sage

    from sage.misc.package import list_packages
    list_packages('optional', local=True)     # Lists optional packages
    list_packages('experimental', local=True) # Lists experimental packages

.. note::

    Use ``local=True`` to avoid network lookups to PyPI, which may fail
    for some Sage package names.

Installation from Source
~~~~~~~~~~~~~~~~~~~~~~~~

When installing from source, you can enable specific optional packages
during the build process. For example, to enable the `bliss` package:

.. code-block:: console

    ./configure --enable-bliss=yes
    make

To see all available configure options:

.. code-block:: console

    ./configure --help

Pip-Installable Optional Packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some optional SageMath packages can be installed using pip:

.. code-block:: console

    ./sage -pip install <package_name>

For example, to install the `graph-genus` package (which provides faster
graph genus algorithms):

.. code-block:: console

    ./sage -pip install graph-genus

This method is particularly useful for Python packages that are not
distributed as standard Sage packages.

Verifying Installation
~~~~~~~~~~~~~~~~~~~~~~

To verify that an optional package was installed correctly, you can check
Sage's package records:

.. code-block:: sage

    from sage.misc.package import list_packages
    pkgs = list_packages('optional', local=True)
    pkgs['bliss'].is_installed()   # True if installed

To check whether the functionality is actually available, use the
corresponding ``sage.features`` check. For example, for `bliss`:

.. code-block:: sage

    from sage.features.bliss import Bliss
    Bliss().is_present()   # True if the feature is available

For instructions on installing optional packages using your system's
package manager (Homebrew, apt, dnf, etc.), see the individual
package pages linked below. Each package page includes the correct
package names for your distribution.

Note
~~~~

Some optional packages may have additional system dependencies. See 
:ref:`All External Packages <spkg>` for more information.

Features
--------

.. toctree::
   :maxdepth: 1

   sage/features
   sage/features/join_feature
   sage/features/all
   sage/features/build_feature
   sage/features/sagemath
   sage/features/pkg_systems
   sage/features/bliss
   sage/features/brial
   sage/features/coxeter3
   sage/features/csdp
   sage/features/databases
   sage/features/dvipng
   sage/features/ffmpeg
   sage/features/four_ti_2
   sage/features/gap
   sage/features/graph_generators
   sage/features/graphviz
   sage/features/imagemagick
   sage/features/interfaces
   sage/features/internet
   sage/features/kenzo
   sage/features/latex
   sage/features/latte
   sage/features/libbraiding
   sage/features/libhomfly
   sage/features/lrs
   sage/features/mcqd
   sage/features/meataxe
   sage/features/mip_backends
   sage/features/normaliz
   sage/features/pandoc
   sage/features/polymake
   sage/features/rankwidth
   sage/features/rubiks
   sage/features/sirocco
   sage/features/singular
   sage/features/tdlib
   sage/features/topcom


Distribution Packages of the Sage Library
-----------------------------------------

.. include:: index_sagemath.rst


Experimental Packages
---------------------

Some packages that provide additional functionality are marked as
"experimental".  Developers are needed in order to improve the
integration of these packages into the Sage distribution.

.. include:: index_experimental.rst


All External Packages
---------------------

.. toctree::
   :maxdepth: 1

   index_alph
