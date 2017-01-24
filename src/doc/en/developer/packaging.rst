.. _chapter-packaging:

==========================
Packaging Third-Party Code
==========================

One of the mottoes of the Sage project is to not reinvent the wheel: If
an algorithm is already implemented in a well-tested library then
consider incorporating that library into Sage. The current list of
available packages are the subdirectories of ``SAGE_ROOT/build/pkgs/``.
The installation of packages is done through a bash script located in
``SAGE_ROOT/build/bin/sage-spkg``. This script is typically invoked by
giving the command::

    [user@localhost]$ sage -i <options> <package name>...

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


.. _section-directory-structure:

Directory Structure
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
   ``FFLAS-FFPACK`` is called ``fflas_ffpack`` in Sage and ``path.py``
   is renamed ``pathpy`` in Sage.

As an example, let us consider a hypothetical FoO project. They
(upstream) distribute a tarball ``FoO-1.3.tar.gz`` (that will be
automatically placed in ``SAGE_ROOT/upstream`` during the installation
process). To package it in Sage, we create a subdirectory containing as
a minimum the following files::

    SAGE_ROOT/build/pkgs/foo
    |-- checksums.ini
    |-- dependencies
    |-- package-version.txt
    |-- spkg-install
    |-- SPKG.txt
    `-- type

The following are some additional files which can be added::

    SAGE_ROOT/build/pkgs/foo
    |-- patches
    |   |-- bar.patch
    |   `-- baz.patch
    |-- spkg-check
    `-- spkg-src

We discuss the individual files in the following sections.


Package type
------------

The file ``type`` should contain a single word, which is either
``standard``, ``optional`` or ``experimental``.
See :ref:`section-package-types` for the meaning of these types.


.. _section-spkg-install:

Build and install scripts
-------------------------

The ``spkg-build`` and ``spkg-install`` files are shell scripts or
Python scripts which build and install the package.
In the best case, the upstream project can simply be installed by the
usual configure / make / make install steps. In that case, the build
script would simply consist of::

    #!/usr/bin/env bash

    cd src

    ./configure --prefix="$SAGE_LOCAL" --libdir="$SAGE_LOCAL/lib"
    if [ $? -ne 0 ]; then
        echo >&2 "Error configuring PACKAGE_NAME."
        exit 1
    fi

    $MAKE
    if [ $? -ne 0 ]; then
        echo >&2 "Error building PACKAGE_NAME."
        exit 1
    fi

The install script would consist of::

    #!/usr/bin/env bash

    cd src
    $MAKE install
    if [ $? -ne 0 ]; then
        echo >&2 "Error installing PACKAGE_NAME."
        exit 1
    fi

Note that the top-level directory inside the tarball is renamed to
``src`` before calling the ``spkg-build`` and ``spkg-install``
scripts, so you can just use ``cd src`` instead of ``cd foo-1.3``.

If there is any meaningful documentation included but not installed by
``make install``, then you can add something like the following to
install it::

    if [ "$SAGE_SPKG_INSTALL_DOCS" = yes ] ; then
        $MAKE doc
        if [ $? -ne 0 ]; then
            echo >&2 "Error building PACKAGE_NAME docs."
            exit 1
        fi
        mkdir -p "$SAGE_SHARE/doc/PACKAGE_NAME"
        cp -R doc/* "$SAGE_SHARE/doc/PACKAGE_NAME"
    fi

Many packages currently do not separate the build and install steps
and only provide a ``spkg-install`` file that does both.  The
separation is useful in particular for root-owned install hierarchies:

- If ``spkg-build`` exists, it is first called, followed by
  ``$SAGE_SUDO spkg-install``.

- Otherwise, only ``spkg-install`` is called (without ``$SAGE_SUDO``).
  Such packages would prefix all commands in ``spkg-install`` that
  write into the installation hierarchy with ``$SAGE_SUDO``.


.. _section-spkg-check:

Self-Tests
----------

The ``spkg-check`` file is an optional, but highly recommended, script
to run self-tests of the package. It is run after building and
installing if the ``SAGE_CHECK`` environment variable is set, see the
Sage installation guide. Ideally, upstream has some sort of tests suite
that can be run with the standard ``make check`` target. In that case,
the ``spkg-check`` script would simply contain::

    #!/usr/bin/env bash

    cd src
    $MAKE check


.. _section-spkg-SPKG-txt:

The SPKG.txt File
-----------------

The ``SPKG.txt`` file should follow this pattern::

     = PACKAGE_NAME =

     == Description ==

     What does the package do?

     == License ==

     What is the license? If non-standard, is it GPLv3+ compatible?

     == Upstream Contact ==

     Provide information for upstream contact.

     == Dependencies ==

     Put a bulleted list of dependencies here:

     * python
     * readline

     == Special Update/Build Instructions ==

     If the tarball was modified by hand and not via a spkg-src
     script, describe what was changed.


with ``PACKAGE_NAME`` replaced by the package name. Legacy
``SPKG.txt`` files have an additional changelog section, but this
information is now kept in the git repository.


.. _section-dependencies:

Package dependencies
--------------------

Many packages depend on other packages. Consider for example the
``eclib`` package for elliptic curves. This package uses the libraries
PARI, NTL and FLINT. So the following is the ``dependencies`` file
for ``eclib``::

    pari ntl flint

    ----------
    All lines of this file are ignored except the first.
    It is copied by SAGE_ROOT/build/make/install into SAGE_ROOT/build/make/Makefile.

If there are no dependencies, you can use ::

    # no dependencies

    ----------
    All lines of this file are ignored except the first.
    It is copied by SAGE_ROOT/build/make/install into SAGE_ROOT/build/make/Makefile.

There are actually two kinds of dependencies: there are normal
dependencies and order-only dependencies, which are weaker. The syntax
for the ``dependencies`` file is ::

    normal dependencies | order-only dependencies

If there is no ``|``, then all dependencies are normal.

- If package A has an **order-only dependency** on B, it simply means
  that B must be built before A can be built. The version of B does not
  matter, only the fact that B is installed matters.
  This should be used if the dependency is purely a build-time
  dependency (for example, a dependency on Python simply because the
  ``spkg-install`` file is written in Python).

- If A has a **normal dependency** on B, it means additionally that A
  should be rebuilt every time that B gets updated. This is the most
  common kind of dependency. A normal dependency is what you need for
  libraries: if we upgrade NTL, we should rebuild everything which
  uses NTL.

In order to check that the dependencies of your package are likely
correct, the following command should work without errors::

    [user@localhost]$ make distclean && make base && make PACKAGE_NAME

Finally, note that standard packages should only depend on standard
packages and optional packages should only depend on standard or
optional packages.


.. _section-spkg-patching:

Patching Sources
----------------

Actual changes to the source code must be via patches, which should be placed
in the ``patches/`` directory, and must have the ``.patch`` extension. GNU
patch is distributed with Sage, so you can rely on it being available. Patches
must include documentation in their header (before the first diff hunk), and
must have only one "prefix" level in the paths (that is, only one path level
above the root of the upstream sources being patched).  So a typical patch file
should look like this::

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
For example, considering the layout::

    SAGE_ROOT/build/pkgs/foo
    |-- patches
    |   |-- solaris
    |   |   |-- solaris.patch 
    |   |-- bar.patch
    |   `-- baz.patch

The patches ``bar.patch`` and ``baz.patch`` are applied to the unpacked
upstream sources in ``src/`` before running ``spkg-install``.  To conditionally
apply the patch for Solaris the ``spkg-install`` should contain a section like
this::

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
  ``Makefile`` because it "hard-codes" some paths or compiler flags::

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
  from the command-line.  Just use the following in ``spkg-install``::

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

- Generate the patches against the ``upstream`` branch::

      rm -Rf SAGE_ROOT/build/pkgs/PACKAGE/patches
      mkdir SAGE_ROOT/build/pkgs/PACKAGE/patches
      git format-patch -o SAGE_ROOT/build/pkgs/PACKAGE/patches/ upstream
  
- Optionally, create an ``spkg-src`` file in the Sage package's
  directory that regenerates the patch directory using the above
  commmands.

- When a new upstream version becomes available, merge (or import) it
  into ``upstream``, then create a new branch and rebase in on top of
  the updated upstream::

      git checkout sage_package_OLDVERSION
      git checkout -b sage_package_NEWVERSION
      git rebase upstream

  Then regenerate the patches.


.. _section-spkg-src:

Modified Tarballs
-----------------

The ``spkg-src`` file is optional and only to document how the upstream
tarball was changed. Ideally it is not modified, then there would be no
``spkg-src`` file present either.

However, if you really must modify the upstream tarball then it is
recommended that you write a script, called ``spkg-src``, that makes the
changes. This not only serves as documentation but also makes it easier
to apply the same modifications to future versions.


.. _section-spkg-versioning:

Package Versioning
------------------

The ``package-version.txt`` file containts just the version. So if
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

Checksums
---------

The ``checksums.ini`` file contains the filename pattern of the
upstream tarball (without the actual version) and its checksums. So if
upstream is ``$SAGE_ROOT/upstream/FoO-1.3.tar.gz``, create a new file
``$SAGE_ROOT/build/pkgs/foo/checksums.ini`` containing only::

    tarball=FoO-VERSION.tar.gz

Sage internally replaces the ``VERSION`` substring with the content of
``package-version.txt``. To recompute the checksums, run::

    [user@localhost]$ sage --package fix-checksum foo

which will modify the ``checksums.ini`` file with the correct
checksums.


Utility script to create package
================================

Assuming that you have downloaded
``$SAGE_ROOT/upstream/FoO-1.3.tar.gz``, you can use::

    [user@localhost]$ sage --package create foo --version 1.3 --tarball FoO-VERSION.tar.gz --type experimental

to create ``$SAGE_ROOT/build/pkgs/foo/package-version.txt``,
``checksums.ini``, and ``type`` in one step.


.. _section-manual-build:

Building the package
====================

At this stage you have a new tarball that is not yet distributed with
Sage (``FoO-1.3.tar.gz`` in the example of section
:ref:`section-directory-structure`). Now you need to manually place it
in the ``SAGE_ROOT/upstream/`` directory and run
``sage --fix-pkg-checksums`` if you have not done that yet.

Now you can install the package using::

    [user@localhost]$ sage -i package_name

or::

    [user@localhost]$ sage -f package_name

to force a reinstallation. If your package contains a ``spkg-check``
script (see :ref:`section-spkg-check`) it can be run with::

    [user@localhost]$ sage -i -c package_name

or::

    [user@localhost]$ sage -f -c package_name

If all went fine, open a ticket, put a link to the original tarball in
the ticket and upload a branch with the code under
``SAGE_ROOT/build/pkgs``.


.. _section-inclusion-procedure:

Inclusion Procedure for New and Updated Packages
================================================

Packages that are not part of Sage will first become optional or
experimental (the latter if they will not build on all supported
systems). After they have been in optional for some time without
problems they can be proposed to be included as standard packages in
Sage.

To propose a package for optional/experimental inclusion please open a
trac ticket with the respective ``Component:`` field set to either
``packages:experimental`` or ``packages:optional``. The associated code
requirements are described in the following sections.

After the ticket was reviewed and included, optional packages stay in
that status for at least a year, after which they can be proposed to be
included as standard packages in Sage. For this a trac ticket is opened
with the ``Component:`` field set to ``packages:standard``. Then make
a proposal in the Google Group ``sage-devel``.

Upgrading packages to new upstream versions or with additional patches
includes opening a ticket in the respective category too, as described
above.

License Information
-------------------

If you are patching a standard Sage spkg, then you should make sure that
the license information for that package is up-to-date, both in its
``SPKG.txt`` file and in the file ``SAGE_ROOT/COPYING.txt``.  For
example, if you are producing an spkg which upgrades the vanilla source
to a new version, check whether the license changed between versions.

Prerequisites for New Standard Packages
---------------------------------------

For a package to become part of Sage's standard distribution, it
must meet the following requirements:

- **License**. For standard packages, the license must be compatible
  with the GNU General Public License, version 3. The Free Software
  Foundation maintains a long list of `licenses and comments about
  them <http://www.gnu.org/licenses/license-list.html>`_.

- **Build Support**. The code must build on all the `fully supported
  platforms
  <http://wiki.sagemath.org/SupportedPlatforms#Fully_supported>`_.

  A standard package should also work on all the platforms where Sage
  is `expected to work
  <http://wiki.sagemath.org/SupportedPlatforms#Expected_to_work>`_ and
  on which Sage `almost works
  <http://wiki.sagemath.org/SupportedPlatforms#Almost_works>`_ but
  since we don't fully support these platforms and often lack the
  resources to test on them, you are not expected to confirm your
  packages works on those platforms.

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
  :ref:`chapter-sage-trac`.


