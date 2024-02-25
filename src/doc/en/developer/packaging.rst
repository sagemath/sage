.. highlight:: shell-session

.. _chapter-packaging:

===================================
Packaging Third-Party Code for Sage
===================================

One of the mottoes of the Sage project is to not reinvent the wheel: If
an algorithm is already implemented in a well-tested library then
consider incorporating that library into Sage. The current list of
available packages are the subdirectories of ``SAGE_ROOT/build/pkgs/``.
The installation of packages is done through a bash script located in
``SAGE_ROOT/build/bin/sage-spkg``. This script is typically invoked by
giving the command::

    [alice@localhost sage]$ sage -i <options> <package name>...

options can be:

- -f: install a package even if the same version is already installed
- -s: do not delete temporary build directory
- -c: after installing, run the test suite for the spkg. This should
  override the settings of ``SAGE_CHECK`` and ``SAGE_CHECK_PACKAGES``.
- -d: only download the package

The section :ref:`section-directory-structure` describes the structure
of each individual package in ``SAGE_ROOT/build/pkgs``. In section
:ref:`section-manual-build` we see how you can install and test a new
spkg that you or someone else wrote. Finally,
:ref:`section-inclusion-procedure` explains how to submit a new package
for inclusion in the Sage source code.


.. _section-package-types:

Package types
=============

Not all packages are built by default, they are divided into standard,
optional and experimental ones:

- **standard** packages are built by default. For a few packages,
  ``configure`` checks whether they are available from the system,
  in which case the build of those packages is skipped.
  Standard packages have stringent quality requirements:
  they should work on all supported platforms. In order
  for a new standard package to be accepted, it should have been
  optional for a while, see :ref:`section-inclusion-procedure`.

- **optional** packages are subject to the same requirements, they
  should also work on all supported platforms. If there are
  :ref:`optional doctests <section-optional-doctest-flag>` in the Sage
  library, those tests must pass.
  Note that optional packages are not tested as much as standard
  packages, so in practice they might break more often than standard
  packages.

- for **experimental** packages, the bar is much lower: even if there are
  some problems, the package can still be accepted.


.. _section-package-source-types:

Package source types
--------------------

Orthogonal to the division by package types, a package has exactly one of
the following source types:

#. A ``normal`` package:

   - comes from the tarball named in the required file ``checksums.ini`` and
     hosted on the Sage mirrors;

   - its version number is defined by the required file ``package-version.txt``;

   - Sage installs the package using build and install scripts
     (see :ref:`section-spkg-install`);

   - Sage records the version number of the package installed using a file in
     ``$SAGE_LOCAL/var/lib/sage/installed/`` and will rerun the installation
     if ``package-version.txt`` changes.

#. A ``wheel`` package:

   - comes from the wheel file named in the required file ``checksums.ini``
     and hosted on the Sage mirrors;

   - per policy, only platform-independent wheels are allowed, i.e.,
     ``*-none-any.whl`` files;

   - its version number is defined by the required file ``package-version.txt``;

   - no build and install scripts are needed
     (with one exception: the package :ref:`spkg_pip` installs itself from
     its wheel using a custom install script);

   - Sage records the version number of the package installed using a file in
     ``$SAGE_LOCAL/var/lib/sage/installed/`` and will rerun the installation
     if ``package-version.txt`` changes.

#. A ``pip`` package:

   - is obtained directly from https://pypi.org/;

   - the version to be installed is determined using the required file
     ``requirements.txt`` -- in its simplest form, this file just
     contains the name of the package (more details at
     https://pip.pypa.io/en/stable/user_guide/#requirements-files);

   - Sage installs the package using the ``pip`` package manager;

   - Sage delegates the recording of installed package version numbers to it;

   - by policy, no ``standard`` package is allowed to be a ``pip`` package.

#. A ``script`` package:

   - is not associated with a tarball;

   - the file ``package-version.txt`` is optional;

   - installing the package runs the installation script ``spkg-install`` or
     ``spkg-install.in`` (see :ref:`section-spkg-install`);

   - Sage records the version number of the package installed using a file in
     ``$SAGE_LOCAL/var/lib/sage/installed/`` and will rerun the installation
     if ``package-version.txt`` changes.

#. A ``dummy`` package:

   - is only used for recording the names of equivalent system packages;

   - there is no ``spkg-install`` script, and attempts to install the package
     using Sage will give an error message.

To summarize: the package source type is determined as follows: if
there is a file ``requirements.txt``, it is a ``pip`` package. If not,
then if there is a ``checksums.ini`` file, it is ``normal`` or ``wheel``.
Otherwise, if it has an ``spkg-install`` or ``spkg-install.in`` script,
it is a ``script`` package, and if it does not, then it is a ``dummy`` package.


.. _section-directory-structure:

Directory structure
===================

Third-party packages in Sage consist of two parts:

#. The tarball as it is distributed by the third party, or as close as
   possible. Valid reasons for modifying the tarball are deleting
   unnecessary files to keep the download size manageable,
   regenerating auto-generated files or changing the directory structure
   if necessary. In certain cases, you may need to (additionally) change
   the filename of the tarball.
   In any case, the actual code must be unmodified: if you need to
   change the sources, add a :ref:`patch <section-spkg-patching>`
   instead. See also :ref:`section-spkg-src` for automating the
   modifications to the upstream tarball.

#. The build scripts and associated files are in a subdirectory
   ``SAGE_ROOT/build/pkgs/<package>``, where you replace ``<package>``
   with a lower-case version of the upstream project name. If the
   project name contains characters which are not alphanumeric
   and are not an underscore, those characters should be removed
   or replaced by an underscore. For example, the project
   ``FFLAS-FFPACK`` is called ``fflas_ffpack`` in Sage.

As an example, let us consider a hypothetical FoO project. They
(upstream) distribute a tarball ``FoO-1.3.tar.gz`` (that will be
automatically placed in ``SAGE_ROOT/upstream`` during the installation
process). To package it in Sage, we create a subdirectory containing as
a minimum the following files:

.. CODE-BLOCK:: text

    SAGE_ROOT/build/pkgs/foo
    |-- checksums.ini
    |-- dependencies
    |-- package-version.txt
    |-- spkg-install.in
    |-- SPKG.rst
    `-- type

The following are some additional files which can be added:

.. CODE-BLOCK:: text

    SAGE_ROOT/build/pkgs/foo
    |-- distros
    |   |-- platform1.txt
    |   `-- platform2.txt
    |-- has_nonfree_dependencies
    |-- huge
    |-- patches
    |   |-- bar.patch
    |   `-- baz.patch
    |-- spkg-check.in
    |-- spkg-configure.m4
    |-- spkg-src
    `-- trees.txt

We discuss the individual files in the following sections.


Package type
------------

The file ``type`` should contain a single word, which is either
``standard``, ``optional`` or ``experimental``.
See :ref:`section-package-types` for the meaning of these types.


.. _section-spkg-install:

Build and install scripts of normal packages
--------------------------------------------

The ``spkg-build.in`` and ``spkg-install.in`` files are templates for
``bash`` scripts ``spkg-build`` and ``spkg-install``, which build
and/or install the package.

The ``*.in`` script templates should *not* be prefixed with a shebang
line (``#!...``) and should not have the executable bit set in their
permissions.  These are added automatically when generating the
scripts, along with some additional boilerplate, when the package is
installed.

The ``spkg-build.in`` and ``spkg-install.in`` files in the Sage source
tree need only focus on the specific steps for building and installing
that package.  If no ``spkg-build.in`` exists, then the
``spkg-install.in`` is responsible for both steps, though separating
them is encouraged where possible.

It is also possible to include similar script templatess named
``spkg-preinst.in`` or ``spkg-postinst.in`` to run additional steps
before or after the package has been installed into
``$SAGE_LOCAL``. It is encouraged to put steps which modify already
installed files in a separate ``spkg-postinst.in`` script template
rather than combining them with ``spkg-install.in``.  This is because
since :trac:`24106`, ``spkg-install`` does not necessarily install
packages directly to ``$SAGE_LOCAL``.  However, by the time
``spkg-postinst`` is run, the installation to ``$SAGE_LOCAL`` is
complete.

In the best case, the upstream project can simply be installed by the
usual configure / make / make install steps. In that case, the
``spkg-build.in`` script template would simply consist of:

.. CODE-BLOCK:: bash

    cd src
    sdh_configure
    sdh_make

See :ref:`section-sdh-helpers` for more on the helper functions
``sdh_configure``, ``sdh_make``, etc.

The ``spkg-install.in`` script template would consist of:

.. CODE-BLOCK:: bash

    cd src
    sdh_make_install

Note that the top-level directory inside the tarball is renamed to
``src`` before calling the ``spkg-build`` and ``spkg-install``
scripts, so you can just use ``cd src`` instead of ``cd foo-1.3``.

If there is any meaningful documentation included but not installed by
``sdh_make_install`` (which calls ``make install``), then you can add
something like the following to install it:

.. CODE-BLOCK:: bash

    if [ "$SAGE_SPKG_INSTALL_DOCS" = yes ] ; then
        sdh_make doc
        sdh_install doc/ "$SAGE_SHARE"/doc/PACKAGE_NAME
    fi

At build time :envvar:`CFLAGS`, :envvar:`CXXFLAGS`, :envvar:`FCFLAGS`,
and :envvar:`F77FLAGS` are usually set to ``-g -O2 -march=native``
(according to `debugging options <../installation/source.html#sage-debug>`_
and whether building
`fat binaries <../installation/source.html#sage-fat-binary>`_).

Slightly modified versions are available:

.. CODE-BLOCK:: bash

    # No ``-march=native``.
    export CFLAGS=$CFLAGS_NON_NATIVE

    # ``-O3`` instead of ``-O2``.
    export CFLAGS=$CFLAGS_O3

    # Use flags as set by the user, possibly empty.
    export CFLAGS=$ORIGINAL_CFLAGS

Likewise for :envvar:`CXXFLAGS`, :envvar:`FCFLAGS`, and :envvar:`F77FLAGS`.

.. note::

    Prior to Sage 9.1, the script templates were called ``spkg-build``,
    ``spkg-install``, etc., without the extension ``.in``.

    Prior to Sage 8.1 the shebang line was included, and the scripts were
    marked executable.  However, this is no longer the case as of
    :trac:`23179`.  Now the scripts in the source tree are deliberately
    written not to be directly executed, and are only made into executable
    scripts when they are copied to the package's build directory.

    Build/install scripts may still be written in Python, but the Python
    code should go in a separate file (e.g. ``spkg-install.py``), and can
    then be executed from the real ``spkg-install.in`` like:

    .. CODE-BLOCK:: text

        exec sage-bootstrap-python spkg-install.py

    or

    .. CODE-BLOCK:: text

        exec python3 spkg-install.py

    In more detail: ``sage-bootstrap-python`` runs a version of Python
    pre-installed on the machine, which is a build prerequisite of Sage.
    Note that ``sage-bootstrap-python`` accepts a wide range of Python
    versions, Python >= 2.6 and >= 3.4, see ``SAGE_ROOT/build/tox.ini``
    for details.  You should only use ``sage-bootstrap-python`` for
    installation tasks that must be able to run before Sage has made
    ``python3`` available.  It must not be used for running ``pip`` or
    ``setup.py`` for any package.

    ``python3`` runs the version of Python managed by Sage (either its
    own installation of Python 3 from an SPKG or a venv over a system
    python3.  You should use this if you are installing a Python package
    to make sure that the libraries are installed in the right place.

    By the way, there is also a script ``sage-python``. This should be
    used at runtime, for example in scripts in ``SAGE_LOCAL/bin`` which
    expect Sage's Python to already be built.

Many packages currently do not separate the build and install steps and only
provide a ``spkg-install.in`` file that does both.  The separation is useful in
particular for root-owned install hierarchies, where something like ``sudo``
must be used to install files.  For this purpose Sage uses an environment
variable ``$SAGE_SUDO``, the value of which may be provided by the developer
at build time,  which should to the appropriate system-specific
``sudo``-like command (if any).  The following rules are then observed:

- If ``spkg-build.in`` exists, the generated script ``spkg-build`` is first
  called, followed by ``$SAGE_SUDO spkg-install``.

- Otherwise, only ``spkg-install`` is called (without ``$SAGE_SUDO``).  Such
  packages should prefix all commands in ``spkg-install.in`` that write into
  the installation hierarchy with ``$SAGE_SUDO``.

Install scripts of script packages
----------------------------------

For script packages, it is also possible to use an install script named ``spkg-install``.
It needs to be an executable shell script; it is not subject to the templating
described in the previous section and will be executed directly from
the build directory.

.. _section-sdh-helpers:

Helper functions
----------------

In the ``spkg-build``, ``spkg-install``, and ``spkg-check`` scripts,
the following functions are available. They are defined in the file
``$SAGE_ROOT/build/bin/sage-dist-helpers``, if you want to look at the
source code.  They should be used to make sure that appropriate
variables are set and to avoid code duplication. These function names
begin with ``sdh_``, which stands for "Sage-distribution helper".

- ``sdh_die MESSAGE``: Exit the build script with the error code of
  the last command if it was non-zero, or with 1 otherwise, and print
  an error message. This is typically used like:

  .. CODE-BLOCK:: bash

       command || sdh_die "Command failed"

  This function can also (if not given any arguments) read the error message
  from stdin. In particular this is useful in conjunction with a heredoc to
  write multi-line error messages:

  .. CODE-BLOCK:: bash

      command || sdh_die << _EOF_
      Command failed.
      Reason given.
      _EOF_

  .. NOTE::

      The other helper functions call ``sdh_die``, so do not use (for
      example) ``sdh_make || sdh_die``: the part of this after
      ``||`` will never be reached.

- ``sdh_check_vars [VARIABLE ...]``: Check that one or more variables
  are defined and non-empty, and exit with an error if any are
  undefined or empty. Variable names should be given without the '$'
  to prevent unwanted expansion.

- ``sdh_configure [...]``: Runs ``./configure`` with arguments
  ``--prefix="$SAGE_LOCAL"``, ``--libdir="$SAGE_LOCAL/lib"``,
  ``--disable-static``, ``--disable-maintainer-mode``, and
  ``--disable-dependency-tracking``. Additional arguments to ``./configure``
  may be given as arguments.

- ``sdh_make [...]``: Runs ``$MAKE`` with the default target.
  Additional arguments to ``$MAKE`` may be given as arguments.

- ``sdh_make_install [...]``: Runs ``$MAKE install`` with DESTDIR
  correctly set to a temporary install directory, for staged
  installations. Additional arguments to ``$MAKE`` may be given as
  arguments. If ``$SAGE_DESTDIR`` is not set then the command is run
  with ``$SAGE_SUDO``, if set.

- ``sdh_setup_bdist_wheel [...]``: Runs ``setup.py bdist_wheel`` with
  the given arguments, as well as additional default arguments used for
  installing packages into Sage.

- ``sdh_pip_install [...]``: The equivalent of running ``pip install``
  with the given arguments, as well as additional default arguments used for
  installing packages into Sage with pip. The last argument must be
  ``.`` to indicate installation from the current directory.

  ``sdh_pip_install`` actually does the installation via ``pip wheel``,
  creating a wheel file in ``dist/``, followed by
  ``sdh_store_and_pip_install_wheel`` (see below).

- ``sdh_pip_editable_install [...]``: The equivalent of running ``pip install -e``
  with the given arguments, as well as additional default arguments used for
  installing packages into Sage with pip. The last argument must be
  ``.`` to indicate installation from the current directory.
  See `pip documentation <https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs>`_
  for more details concerning editable installs.

- ``sdh_pip_uninstall [...]``: Runs ``pip uninstall`` with the given arguments.
  If unsuccessful, it displays a warning.

- ``sdh_store_and_pip_install_wheel .``: The current directory,
  indicated by the required argument ``.``, must have a subdirectory
  ``dist`` containing a unique wheel file (``*.whl``).

  This command (1) moves this wheel file to the
  directory ``$SAGE_SPKG_WHEELS`` (``$SAGE_LOCAL/var/lib/sage/wheels``)
  and then (2) installs the wheel in ``$SAGE_LOCAL``.

  Both of these steps, instead of writing directly into ``$SAGE_LOCAL``,
  use the staging directory ``$SAGE_DESTDIR`` if set; otherwise, they
  use ``$SAGE_SUDO`` (if set).

- ``sdh_install [-T] SRC [SRC...] DEST``: Copies one or more files or
  directories given as ``SRC`` (recursively in the case of
  directories) into the destination directory ``DEST``, while
  ensuring that ``DEST`` and all its parent directories exist.
  ``DEST`` should be a path under ``$SAGE_LOCAL``, generally. For
  ``DESTDIR`` installs, the ``$SAGE_DESTDIR`` path is automatically
  prepended to the destination.

  The ``-T`` option treats ``DEST`` as a normal file instead
  (e.g. for copying a file to a different filename). All directory
  components are still created in this case.

The following is automatically added to each install script, so you
should not need to add it yourself.

- ``sdh_guard``: Wrapper for ``sdh_check_vars`` that checks some
  common variables without which many/most packages won't build
  correctly (``SAGE_ROOT``, ``SAGE_LOCAL``, ``SAGE_SHARE``). This is
  important to prevent installation to unintended locations.

The following are also available, but rarely used.

- ``sdh_cmake [...]``: Runs ``cmake`` in the current directory with
  the given arguments, as well as additional arguments passed to
  cmake (assuming packages are using the GNUInstallDirs module) so
  that ``CMAKE_INSTALL_PREFIX`` and ``CMAKE_INSTALL_LIBDIR`` are set
  correctly.

- ``sdh_preload_lib EXECUTABLE SONAME``: (Linux only -- no-op on other
  platforms.)  Check shared libraries loaded by ``EXECUTABLE`` (may be a
  program or another library) for a library starting with ``SONAME``, and
  if found appends ``SONAME`` to the ``LD_PRELOAD`` environment variable.
  See :trac:`24885`.


.. _spkg-configure.m4:

Allowing for the use of system packages
---------------------------------------

For a number of Sage packages, an already installed system version can
be used instead, and Sage's top-level ``./configure`` script
determines when this is possible. To enable this, a package needs to
have a script called ``spkg-configure.m4``, which can, for example,
determines whether the installed software is recent enough (and
sometimes not too recent) to be usable by Sage. This script is
processed by the `GNU M4 macro processor
<https://www.gnu.org/savannah-checkouts/gnu/m4/manual/m4-1.4.18/m4.html>`_.

Also, if the software for a Sage package is provided by a system
package, the ``./configure`` script can provide that information. To
do this, there must be a directory ``build/pkgs/PACKAGE/distros``
containing files with names like ::

    arch.txt
    conda.txt
    debian.txt
    fedora.txt
    homebrew.txt
    ...

corresponding to different packaging systems. Each system package
should appear on a separate line. If the shell-style variable reference
``${PYTHON_MINOR}`` appears, it is replaced by the minor version of
Python, e.g., 12 for Python 3.12.x. Everything on a line after a ``#``
character is ignored, so comments can be included in the files.

For example, if ``./configure`` detects that the Homebrew packaging
system is in use, and if the current package can be provided by a
Homebrew package called "foo", then the file
``build/pkgs/PACKAGE/distros/homebrew.txt`` should contain the single
line "foo". If ``foo`` is currently uninstalled, then ``./configure``
will print a message suggesting that the user should run ``brew install
foo``. See :ref:`section-equiv-distro-packages` for more on this.

.. IMPORTANT::

    All new standard packages should, when possible, include a
    ``spkg-configure.m4`` script and a populated ``distros``
    directory. There are many examples in ``build/pkgs``, including
    ``build/pkgs/python3`` and ``build/pkgs/suitesparse``, to name a few.

Note that this may not be possible (as of this writing) for some
packages, for example packages installed via pip for use while running
Sage, like ``matplotlib`` or ``scipy``. If a package is installed via
pip for use in a separate process, like ``tox``, then this should be
possible.



.. _section-spkg-check:

Self-tests
----------

The ``spkg-check.in`` file is an optional, but highly recommended,
script template to run self-tests of the package.  The format for the
``spkg-check`` is the same as ``spkg-build`` and ``spkg-install``.  It
is run after building and installing if the ``SAGE_CHECK`` environment
variable is set, see the Sage installation guide. Ideally, upstream
has some sort of test suite that can be run with the standard ``make
check`` target. In that case, the ``spkg-check.in`` script template
would simply contain:

.. CODE-BLOCK:: bash

    cd src
    $MAKE check


.. _section-python:

Python-based packages
---------------------

Python-based packages should declare ``$(PYTHON)`` as a dependency,
and most Python-based packages will also have ``$(PYTHON_TOOLCHAIN)`` as
an order-only dependency, which will ensure that fundamental packages such
as ``pip`` and ``setuptools`` are available at the time of building the package.

The best way to install a Python-based package is to use ``pip``, in which
case the ``spkg-install.in`` script template might just consist of

.. CODE-BLOCK:: bash

    cd src && sdh_pip_install .

Where ``sdh_pip_install`` is a function provided by ``sage-dist-helpers`` that
points to the correct ``pip`` for the Python used by Sage, and includes some
default flags needed for correct installation into Sage.

If ``pip`` will not work for a package but a command like ``python3 setup.py install``
will, you may use ``sdh_setup_bdist_wheel``, followed by
``sdh_store_and_pip_install_wheel .``.

For ``spkg-check.in`` script templates, use ``python3`` rather
than just ``python``.  The paths are set by the Sage build system
so that this runs the correct version of Python.
For example, the ``scipy`` ``spkg-check.in`` file contains the line

.. CODE-BLOCK:: bash

    exec python3 spkg-check.py

All normal Python packages and all wheel packages must have a file ``install-requires.txt``.
If a Python package is available on PyPI, this file must contain the
name of the package as it is known to PyPI.  Optionally,
``install-requires.txt`` can encode version constraints (such as lower
and upper bounds).  The constraints are in the format of the
``install_requires`` key of `setup.cfg
<https://setuptools.readthedocs.io/en/latest/userguide/declarative_config.html>`_
or `setup.py
<https://packaging.python.org/discussions/install-requires-vs-requirements/#id5>`_.

It is strongly recommended to include comments (starting with ``#``)
in the file that explain why a particular lower or upper bound is
warranted or why we wish to include or reject certain versions.

For example:

.. CODE-BLOCK:: bash

    $ cat build/pkgs/sphinx/package-version.txt
    3.1.2.p0
    $ cat build/pkgs/sphinx/install-requires.txt
    # gentoo uses 3.2.1
    sphinx >=3, <3.3

The comments may include links to GitHub Issues/PRs, as in the following example:

.. CODE-BLOCK:: bash

    $ cat build/pkgs/packaging/install-requires.txt
    packaging >=18.0
    # Issue #30975: packaging 20.5 is known to work
    # but we have to silence "DeprecationWarning: Creating a LegacyVersion"

The currently encoded version constraints are merely a starting point.
Developers and downstream packagers are invited to refine the version
constraints based on their experience and tests.  When a package
update is made in order to pick up a critical bug fix from a newer
version, then the lower bound should be adjusted.
Setting upper bounds to guard against incompatible future changes is
a complex topic; see :trac:`33520`.


.. _section-spkg-SPKG-txt:

The SPKG.rst file
-----------------

The ``SPKG.rst`` file should follow this pattern:

.. CODE-BLOCK:: text

     PACKAGE_NAME: One line description
     ==================================

     Description
     -----------

     What does the package do?

     License
     -------

     What is the license? If non-standard, is it GPLv3+ compatible?

     Upstream Contact
     ----------------

     Provide information for upstream contact.  Usually just an URL.

     Dependencies
     ------------

     Only put special dependencies here that are not captured by the
     ``dependencies`` file. Otherwise omit this section.

     Special Update/Build Instructions
     ---------------------------------

     If the tarball was modified by hand and not via an ``spkg-src``
     script, describe what was changed. Otherwise omit this section.


with ``PACKAGE_NAME`` replaced by the SPKG name (= the directory name in
``build/pkgs``).

Legacy ``SPKG.txt`` files have an additional changelog section, but this
information is now kept in the git repository.

.. _section-dependencies:

Package dependencies
--------------------

Many packages depend on other packages. Consider for example the
``eclib`` package for elliptic curves. This package uses the libraries
PARI, NTL and FLINT. So the following is the ``dependencies`` file
for ``eclib``:

.. CODE-BLOCK:: text

    pari ntl flint

    ----------
    All lines of this file are ignored except the first.

For Python packages, common dependencies include ``pip``,
``setuptools``, and ``future``. If your package depends on any of
these, use ``$(PYTHON_TOOLCHAIN)`` instead. For example, here is the
``dependencies`` file for ``configparser``:

.. CODE-BLOCK:: text

    $(PYTHON) | $(PYTHON_TOOLCHAIN)

(See below for the meaning of the ``|``.)

If there are no dependencies, you can use

.. CODE-BLOCK:: text

    # no dependencies

    ----------
    All lines of this file are ignored except the first.

There are actually two kinds of dependencies: there are normal
dependencies and order-only dependencies, which are weaker. The syntax
for the ``dependencies`` file is

.. CODE-BLOCK:: text

    normal dependencies | order-only dependencies

If there is no ``|``, then all dependencies are normal.

- If package A has an **order-only dependency** on B, it simply means
  that B must be built before A can be built. The version of B does not
  matter, only the fact that B is installed matters.
  This should be used if the dependency is purely a build-time
  dependency (for example, a dependency on pip simply because the
  ``spkg-install`` file uses pip).

  Alternatively, you can put the order-only dependencies in a separate
  file ``dependencies_order_only``.

- If A has a **normal dependency** on B, it means additionally that A
  should be rebuilt every time that B gets updated. This is the most
  common kind of dependency. A normal dependency is what you need for
  libraries: if we upgrade NTL, we should rebuild everything which
  uses NTL.

Some packages are only needed for self-tests of a package (``spkg-check``).
These dependencies should be declared in a separate file ``dependencies_check``.

Some dependencies are optional in the sense that they are only
a dependency if they are configured to be installed. These dependencies
should be declared in a separate file ``dependencies_optional``.

In order to check that the dependencies of your package are likely
correct, the following command should work without errors::

    [alice@localhost sage]$ make distclean && make base && make PACKAGE_NAME

Finally, note that standard packages should only depend on standard
packages and optional packages should only depend on standard or
optional packages.


.. _section-spkg-tags:

Package tags
------------

You can mark a package as "huge" by placing an empty file named
``huge`` in the package directory.  For example, the package
``polytopes_db_4d`` is a large database whose compressed tarball has a
size of 9 GB.

For some other packages, we have placed an empty file named
``has_nonfree_dependencies`` in the package directory. This is to
indicate that Sage with this package installed cannot be
redistributed, and also that the package can only be installed after
installing some other, non-free package.

We use these tags in our continuous integration scripts to filter
out packages that we cannot or should not test automatically.


.. _section-trees:

Where packages are installed
----------------------------

The Sage distribution has the notion of several installation trees.

- ``$SAGE_VENV`` is the default installation tree for all Python packages, i.e.,
  normal packages with an ``install-requires.txt``, wheel packages, and pip packages
  with a ``requirements.txt``.

- ``$SAGE_LOCAL`` is the default installation tree for all non-Python packages.

- ``$SAGE_DOCS`` (only set at build time) is an installation tree for the
  HTML and PDF documentation.

By placing a file ``trees.txt`` in the package directory, the installation tree
can be overridden.  For example, ``build/pkgs/python3/trees.txt`` contains the
word ``SAGE_VENV``, and ``build/pkgs/sagemath_doc_html/trees.txt`` contains the
word ``SAGE_DOCS``.


.. _section-spkg-patching:

Patching sources
----------------

Actual changes to the source code must be via patches, which should be placed
in the ``patches/`` directory, and must have the ``.patch`` extension. GNU
patch is distributed with Sage, so you can rely on it being available. Patches
must include documentation in their header (before the first diff hunk), and
must have only one "prefix" level in the paths (that is, only one path level
above the root of the upstream sources being patched).  So a typical patch file
should look like this:

.. CODE-BLOCK:: diff

    Add autodoc_builtin_argspec config option

    Following the title line you can add a multi-line description of
    what the patch does, where you got it from if you did not write it
    yourself, if they are platform specific, if they should be pushed
    upstream, etc...

    diff -dru Sphinx-1.2.2/sphinx/ext/autodoc.py.orig Sphinx-1.2.2/sphinx/ext/autodoc.py
    --- Sphinx-1.2.2/sphinx/ext/autodoc.py.orig  2014-03-02 20:38:09.000000000 +1300
    +++ Sphinx-1.2.2/sphinx/ext/autodoc.py  2014-10-19 23:02:09.000000000 +1300
    @@ -1452,6 +1462,7 @@

         app.add_config_value('autoclass_content', 'class', True)
         app.add_config_value('autodoc_member_order', 'alphabetic', True)
    +    app.add_config_value('autodoc_builtin_argspec', None, True)
         app.add_config_value('autodoc_default_flags', [], True)
         app.add_config_value('autodoc_docstring_signature', True, True)
         app.add_event('autodoc-process-docstring')

Patches directly under the ``patches/`` directly are applied automatically
before running the ``spkg-install`` script (so long as they have the ``.patch``
extension).  If you need to apply patches conditionally (such as only on
a specifically platform), you can place those patches in a subdirectory of
``patches/`` and apply them manually using the ``sage-apply-patches`` script.
For example, considering the layout:

.. CODE-BLOCK:: text

    SAGE_ROOT/build/pkgs/foo
    |-- patches
    |   |-- solaris
    |   |   |-- solaris.patch
    |   |-- bar.patch
    |   `-- baz.patch

The patches ``bar.patch`` and ``baz.patch`` are applied to the unpacked
upstream sources in ``src/`` before running ``spkg-install``.  To conditionally
apply the patch for Solaris the ``spkg-install`` should contain a section like
this:

.. CODE-BLOCK:: bash

    if [ $UNAME == "SunOS" ]; then
        sage-apply-patches -d solaris
    fi

where the ``-d`` flag applies all patches in the ``solaris/`` subdirectory of
the main ``patches/`` directory.


.. _section-spkg-patch-or-repackage:

When to patch, when to repackage, when to autoconfiscate
--------------------------------------------------------

- Use unpatched original upstream tarball when possible.

  Sometimes it may seem as if you need to patch a (hand-written)
  ``Makefile`` because it "hard-codes" some paths or compiler flags:

  .. CODE-BLOCK:: diff

      --- a/Makefile
      +++ b/Makefile
      @@ -77,7 +77,7 @@
       # This is a Makefile.
       # Handwritten.

      -DESTDIR = /usr/local
      +DESTDIR = $(SAGE_ROOT)/local
       BINDIR   = $(DESTDIR)/bin
       INCDIR   = $(DESTDIR)/include
       LIBDIR   = $(DESTDIR)/lib

  Don't use patching for that.  Makefile variables can be overridden
  from the command-line.  Just use the following in ``spkg-install``:

  .. CODE-BLOCK:: bash

      $(MAKE) DESTDIR="$SAGE_ROOT/local"

- Check if Debian or another distribution already provides patches
  for upstream.  Use them, don't reinvent the wheel.

- If the upstream Makefile does not build shared libraries,
  don't bother trying to patch it.

  Autoconfiscate the package instead and use the standard facilities
  of Automake and Libtool.  This ensures that the shared library build
  is portable between Linux and macOS.

- If you have to make changes to ``configure.ac`` or other source
  files of the autotools build system (or if you are autoconfiscating
  the package), then you can't use patching; make a :ref:`modified
  tarball <section-spkg-src>` instead.

- If the patch would be huge, don't use patching.  Make a
  :ref:`modified tarball <section-spkg-src>` instead.

- Otherwise, :ref:`maintain a set of patches
  <section-spkg-patch-maintenance>`.


.. _section-spkg-patch-maintenance:

How to maintain a set of patches
--------------------------------

We recommend the following workflow for maintaining a set of patches.

- Fork the package and put it on a public git repository.

  If upstream has a public version control repository, import it from
  there.  If upstream does not have a public version control
  repository, import the current sources from the upstream tarball.
  Let's call the branch ``upstream``.

- Create a branch for the changes necessary for Sage, let's call it
  ``sage_package_VERSION``, where ``version`` is the upstream version
  number.

- Make the changes and commit them to the branch.

- Generate the patches against the ``upstream`` branch:

  .. CODE-BLOCK:: bash

      rm -Rf SAGE_ROOT/build/pkgs/PACKAGE/patches
      mkdir SAGE_ROOT/build/pkgs/PACKAGE/patches
      git format-patch -o SAGE_ROOT/build/pkgs/PACKAGE/patches/ upstream

- Optionally, create an ``spkg-src`` file in the Sage package's
  directory that regenerates the patch directory using the above
  commands.

- When a new upstream version becomes available, merge (or import) it
  into ``upstream``, then create a new branch and rebase it on top of
  the updated upstream:

  .. CODE-BLOCK:: bash

      git checkout sage_package_OLDVERSION
      git checkout -b sage_package_NEWVERSION
      git rebase upstream

  Then regenerate the patches.


.. _section-spkg-src:

Modified tarballs
-----------------

The ``spkg-src`` file is optional and only to document how the upstream
tarball was changed. Ideally it is not modified, then there would be no
``spkg-src`` file present either.

However, if you really must modify the upstream tarball then it is
recommended that you write a script, called ``spkg-src``, that makes the
changes. This not only serves as documentation but also makes it easier
to apply the same modifications to future versions.


.. _section-spkg-versioning:

Package versioning
------------------

The ``package-version.txt`` file contains just the version. So if
upstream is ``FoO-1.3.tar.gz`` then the package version file would only
contain ``1.3``.

If the upstream package is taken from some revision other than a stable
version or if upstream doesn't have a version number, you should use the
date at which the revision is made. For example, the
``database_stein_watkins`` package with version ``20110713`` contains
the database as of 2011-07-13. Note that the date should refer to the
*contents* of the tarball, not to the day it was packaged for Sage.
This particular Sage package for ``database_stein_watkins`` was created
in 2014, but the data it contains was last updated in 2011.

If you apply any patches, or if you made changes to the upstream tarball
(see :ref:`section-directory-structure` for allowable changes),
then you should append a ``.p0`` to the version to indicate that it's
not a vanilla package.

Additionally, whenever you make changes to a package *without* changing
the upstream tarball (for example, you add an additional patch or you
fix something in the ``spkg-install`` file), you should also add or
increase the patch level. So the different versions would
be ``1.3``, ``1.3.p0``, ``1.3.p1``, ...
The change in version number or patch level will trigger
re-installation of the package, such that the changes are taken into
account.


.. _section-spkg-checksums:

Checksums and tarball names
---------------------------

The ``checksums.ini`` file contains the filename pattern of the
upstream tarball (without the actual version) and its checksums. So if
upstream is ``$SAGE_ROOT/upstream/FoO-1.3.tar.gz``, create a new file
``$SAGE_ROOT/build/pkgs/foo/checksums.ini`` containing only:

.. CODE-BLOCK:: bash

    tarball=FoO-VERSION.tar.gz

Sage internally replaces the ``VERSION`` substring with the content of
``package-version.txt``.


.. _section-spkg-upstream-urls:

Upstream URLs
-------------

In addition to these fields in ``checksums.ini``, the optional field
``upstream_url`` holds an URL to the upstream package archive.

The Release Manager uses the information in ``upstream_url`` to
download the upstream package archive and to make it available on the
Sage mirrors when a new release is prepared.  On GitHub PRs
upgrading a package, the PR description should no longer contain
the upstream URL to avoid duplication of information.

Note that, like the ``tarball`` field, the ``upstream_url`` is a
template; the substring ``VERSION`` is substituted with the actual
version. It can also be written as ``${VERSION}``, and it is possible
to refer to the dot-separated components of a version by ``VERSION_MAJOR``,
``VERSION_MINOR``, and ``VERSION_MICRO``.

For Python packages available from PyPI, you should use an
``upstream_url`` from ``pypi.io``, which follows the format

.. CODE-BLOCK:: bash

    upstream_url=https://pypi.io/packages/source/m/matplotlib/matplotlib-VERSION.tar.gz

Developers who wish to test a package update from a PR branch before
the archive is available on a Sage mirror. Sage falls back to
downloading package tarballs from the ``upstream_url`` after trying all
Sage mirrors. (This can be disabled by using ``./configure
--disable-download-from-upstream-url``.)  To speed up this process,
you can trim ``upstream/mirror_list`` to fewer mirrors.


.. _section-sage-package-command:

Utility script to create and maintain packages
==============================================

The command ``sage --package`` offers a range of functionality for
creating and maintaining packages of the Sage distribution.


Creating packages
-----------------

Assuming that you have downloaded
``$SAGE_ROOT/upstream/FoO-1.3.tar.gz``, you can use::

    [alice@localhost sage]$ sage --package create foo                     \
                                             --version 1.3                \
                                             --tarball FoO-VERSION.tar.gz \
                                             --type experimental

to create ``$SAGE_ROOT/build/pkgs/foo/package-version.txt``,
``checksums.ini``, and ``type`` in one step.

You can skip the manual downloading of the upstream tarball by using
the additional argument ``--upstream-url``.  This command will also
set the ``upstream_url`` field in ``checksums.ini`` described above.

For Python packages available from PyPI, you can use::

    [alice@localhost sage]$ sage --package create scikit_spatial --pypi   \
                                             --type optional

This automatically downloads the most recent version from PyPI and also
obtains most of the necessary information by querying PyPI.

The ``dependencies`` file may need editing (watch out for warnings regarding
``--no-deps`` that Sage issues during installation of the package!).
Also you may want to set lower and upper bounds for acceptable package versions
in the file ``install-requires.txt``.

By default, when the package is available as a platform-independent
wheel, the ``sage --package`` creates a wheel package. To create a normal package
instead (for example, when the package requires patching), you can use::

    [alice@localhost sage]$ sage --package create scikit_spatial --pypi   \
                                             --source normal              \
                                             --type optional

To create a pip package rather than a normal or wheel package, you can use::

    [alice@localhost sage]$ sage --package create scikit_spatial --pypi   \
                                             --source pip                 \
                                             --type optional

When the package already exists, ``sage --package create`` overwrites it.


Updating packages to a new version
----------------------------------

A package that has the ``upstream_url`` information can be updated by
simply typing::

    [alice@localhost sage]$ sage --package update numpy 3.14.59

which will automatically download the archive and update the
information in ``build/pkgs/numpy/``.

For Python packages available from PyPI, there is another shortcut::

    [alice@localhost sage]$ sage --package update-latest matplotlib
    Updating matplotlib: 3.3.0 -> 3.3.1
    Downloading tarball to ...matplotlib-3.3.1.tar.bz2
    [...............................................................]

If you pass the switch ``--commit``, the script will run ``git commit``
for you.

If you prefer to make update a package ``foo`` by making manual
changes to the files in ``build/pkgs/foo``, you will need to run::

    [alice@localhost sage]$ sage --package fix-checksum foo

which will modify the ``checksums.ini`` file with the correct
checksums.


Obtaining package metrics
-------------------------

The command ``sage --package metrics`` computes machine-readable
aggregated metrics for all packages in the Sage distribution or a
given list of packages::

    [alice@localhost sage]$ sage --package metrics
    has_file_distros_arch_txt=181
    has_file_distros_conda_txt=289
    has_file_distros_debian_txt=172
    has_file_distros_fedora_txt=183
    has_file_distros_gentoo_txt=211
    has_file_distros_homebrew_txt=95
    has_file_distros_macports_txt=173
    has_file_distros_nix_txt=72
    has_file_distros_opensuse_txt=206
    has_file_distros_slackware_txt=32
    has_file_distros_void_txt=221
    has_file_patches=63
    has_file_spkg_check=106
    has_file_spkg_configure_m4=262
    has_file_spkg_install=322
    has_tarball_upstream_url=291
    line_count_file_patches=31904
    line_count_file_spkg_check=585
    line_count_file_spkg_configure_m4=3337
    line_count_file_spkg_install=4342
    packages=442
    type_base=1
    type_experimental=18
    type_optional=151
    type_standard=272

Developers can use these metrics to monitor the complexity and quality
of the Sage distribution. Here are some examples:

- ``has_file_patches`` indicates how many packages have non-empty
  ``patches/`` directories, and ``line_count_file_patches`` gives
  the total number of lines in the patch files.

  Ideally, we would not have to carry patches for a
  package. For example, updating patches when a new upstream version
  is released can be a maintenance burden.

  Developers can help by working with the upstream maintainers of the
  package to prepare a new version that requires fewer or smaller
  patches, or none at all.

- ``line_count_spkg_install`` gives the total number of lines in
  ``spkg-install`` or ``spkg-install.in`` files; see
  :ref:`section-spkg-install`.

  When we carry complex ``spkg-install.in`` scripts for normal packages,
  it may indicate that the upstream package's build and installation
  scripts should be improved.

  Developers can help by working with the upstream maintainers of the
  package to prepare an improved version.

- ``has_file_spkg_check`` indicates how many packages have an
  ``spkg-check`` or ``spkg-check.in`` file; see :ref:`section-spkg-check`.

- ``has_file_spkg_configure_m4`` indicates how many packages
  are prepared to check for an equivalent system package, and
  ``has_file_distros_arch_txt``, ``has_file_distros_conda_txt``
  etc. count how many packages provide the corresponding system package
  information.


.. _section-manual-build:

Building the package
====================

At this stage you have a new tarball that is not yet distributed with
Sage (``FoO-1.3.tar.gz`` in the example of section
:ref:`section-directory-structure`).

Now you can install the package using::

    [alice@localhost sage]$ sage -i package_name

or::

    [alice@localhost sage]$ sage -f package_name

to force a reinstallation. If your package contains a ``spkg-check``
script (see :ref:`section-spkg-check`) it can be run with::

    [alice@localhost sage]$ sage -i -c package_name

or::

    [alice@localhost sage]$ sage -f -c package_name

If all went fine, open a PR with the code under
``SAGE_ROOT/build/pkgs``.


.. _section-inclusion-procedure:

Inclusion procedure for new and updated packages
================================================

Packages that are not part of Sage will first become optional or
experimental (the latter if they will not build on all supported
systems). After they have been in optional for some time without
problems they can be proposed to be included as standard packages in
Sage.

To propose a package for optional/experimental inclusion please open a GitHub
PR added with labels ``c: packages: experimental`` or ``c: packages:
optional``. The associated code requirements are described in the following
sections.

After the PR was reviewed and included, optional packages stay in
that status for at least a year, after which they can be proposed to be
included as standard packages in Sage. For this a GitHub PR is opened
with the label ``c: packages: standard``. Then make
a proposal in the Google Group ``sage-devel``.

Upgrading packages to new upstream versions or with additional patches
includes opening a PR in the respective category too, as described
above.

License information
-------------------

If you are patching a standard Sage spkg, then you should make sure that
the license information for that package is up-to-date, both in its
``SPKG.rst`` or ``SPKG.txt`` file and in the file ``SAGE_ROOT/COPYING.txt``.  For
example, if you are producing an spkg which upgrades the vanilla source
to a new version, check whether the license changed between versions.

If an upstream tarball of a package cannot be redistributed for license
reasons, rename it to include the string ``do-not-distribute``.  This
will keep the release management scripts from uploading it to the Sage mirrors.

Sometimes an upstream tarball contains some distributable parts using
a free software license and some non-free parts.  In this case, it can
be a good solution to make a custom tarball consisting of only the free
parts; see :ref:`section-spkg-src` and the ``giac`` package as an example.


Prerequisites for new standard packages
---------------------------------------

For a package to become part of Sage's standard distribution, it
must meet the following requirements:

- **License**. For standard packages, the license must be compatible
  with the GNU General Public License, version 3. The Free Software
  Foundation maintains a long list of `licenses and comments about
  them <http://www.gnu.org/licenses/license-list.html>`_.

- **Build Support**. The code must build on all the fully supported
  platforms (Linux, macOS); see :ref:`chapter-portability_testing`.
  It must be installed either from source as a normal package,
  or as a Python (platform-independent) wheel package, see
  :ref:`section-package-source-types`.

- **Quality**. The code should be "better" than any other available
  code (that passes the two above criteria), and the authors need to
  justify this. The comparison should be made to both Python and other
  software. Criteria in passing the quality test include:

  - Speed

  - Documentation

  - Usability

  - Absence of memory leaks

  - Maintainable

  - Portability

  - Reasonable build time, size, dependencies

- **Previously an optional package**. A new standard package must have
  spent some time as an optional package. Or have a good reason why
  this is not possible.

- **Refereeing**. The code must be refereed, as discussed in
  :ref:`chapter-github`.


