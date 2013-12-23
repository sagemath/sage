.. _chapter-old-spkg:

=========================
Packaging Old-Style SPKGs
=========================

This chapter explains old-style spkgs; It applies only to legacy
optional spkgs and experimental spkgs. See :ref:`chapter-packaging`
for the modern way of packaging third-party software.


Creating an Old-Style SPKG
==========================

If you are producing code to add new functionality to Sage, you might
consider turning it into a package (an "spkg") instead of a patch
file. If your code is very large (for instance) and should be offered
as an optional download, a package is the right choice. Similarly, if
your code depends on some other optional component of Sage, you should
produce a package. When in doubt, ask for advice on the ``sage-devel``
mailing list.

The abbreviation "spkg" stands for "Sage package". The directory
``SAGE_ROOT/spkg/standard`` contains spkg's. In a source install,
these are all Sage spkg files (actually ``.tar`` or ``.tar.bz2``
files), which are the source code that defines Sage. In a binary
install, some of these may be small placeholder files to save space.

Sage packages are distributed as ``.spkg`` files, which are
``.tar.bz2`` files (or ``tar`` files) but have the extension ``.spkg``
to discourage confusion. Although Sage packages are packed using tar
and/or bzip2, note that ``.spkg`` files contain control information
(installation scripts and metadata) that are necessary for building
and installing them.  When you compile Sage from a source distribution
(or when you run ``sage -i <pkg>`` or ``sage -f <pkg>``), the file
``SAGE_ROOT/src/bin/sage-spkg`` takes care of the unpacking,
compilation, and installation of Sage packages for you. You can type::

    tar -jxvf mypackage-version.spkg

to extract an spkg and see what is inside.  If you want to create a
new Sage package, it is recommended that you start by examining some
existing spkg's. The URL
http://www.sagemath.org/download-packages.html lists spkg's available
for download.


Naming Your SPKG
----------------

Each Sage spkg has a name of the following form::

   BASENAME-VERSION.spkg

``BASENAME`` is the name of the package; it may contain lower-case
letters, numbers, and underscores, but no hyphens.  ``VERSION`` is the
version number; it should start with a number and may contain numbers,
letters, dots, and hyphens; it may end in a string of the form "pNUM",
where "NUM" is a non-negative integer.  If your spkg is a "vanilla"
(unmodified) version of some piece of software, say version 5.3 of
"my-python-package", then ``BASENAME`` would be "my_python_package" --
note the change from hyphens to underscores, because ``BASENAME``
should not contain any hyphens -- and ``VERSION`` would be "5.3".  If
you need to modify the software to use it with Sage (as described
below and in the chapter :ref:`section-old-spkg-patching-overview`),
then ``VERSION`` would be "5.3.p0", the "p0" indicating a patch-level
of 0.  If someone adds more patches, later, this would become "p1",
then "p2", etc.

The string ``VERSION`` must be present.  If you are using a piece
software with no obvious version number, use a date. To give your spkg
a name like this, create a directory called ``BASENAME-VERSION`` and
put your files in that directory -- the next section describes the
directory structure.


Directory Structure
-------------------

Put your files in a directory with a name like ``mypackage-0.1``, as
described above.  If you are porting
another software package, then the directory should contain a
subdirectory ``src/``, containing an unaltered copy of the package.
Every file not in ``src/`` should be under version control, i.e. checked
into an hg repository.

More precisely, the directory should contain the following:

- ``src/``: this directory contains vanilla upstream code, with a few
  exceptions, e.g. when the spkg shipped with Sage is in effect
  upstream, and development on that code base is happening in close
  coordination with Sage.  See John Cremona's  eclib spkg, for
  instance. The directory ``src/`` must not be under revision control.

- ``.hg``, ``.hgignore``, and ``.hgtags``: Old-style spkgs use
  Mercurial for its revision control system. The hidden directory
  ``.hg`` is part of the standard Sage spkg layout.  It contains the
  Mercurial repository for all files not in the ``src/`` directory.
  To create this Mercurial repository from scratch, you should do::

      hg init

  The files ``.hgignore`` and ``.hgtags`` also belong to the Mercurial
  repository.  The file ``.hgtags`` is optional, and is frequently
  omitted.  You should make sure that the file ``.hgignore`` contains
  "src/", since we are not tracking its content.  Indeed, frequently
  this file contains only a single line::

      src/

- ``spkg-install``: this file contains the install script.  See
  :ref:`section-old-spkg-install` for more information and a template.

- ``SPKG.txt``: this file describes the spkg in wiki format.  Each new
  revision needs an updated changelog entry or the spkg will get an
  automatic "needs work" at review time.  See
  :ref:`section-old-spkg-SPKG-txt` for a template.

- ``spkg-check``: this file runs the test suite.  This is somewhat
  optional since not all spkg's have test suites. If possible, do
  create such a script since it helps isolate bugs in upstream
  packages.

- ``patches/``: this directory contains patches to source files in
  ``src/``.  See :ref:`section-old-spkg-patching-overview`.  Patches
  to files in ``src/`` should be applied in ``spkg-install``, and all
  patches must be documented in ``SPKG.txt``, i.e. what they do, if
  they are platform specific, if they should be pushed upstream,
  etc. To ensure that all patched versions of upstream source files
  under ``src/`` are under revision control, the whole directory
  ``patches/`` must be under revision control.

**Never** apply patches to upstream source files under ``src/`` and
then package up an spkg. Such a mixture of upstream source with Sage
specific patched versions is a recipe for confusion. There must be a
**clean separation** between the source provided by the upstream
project and the patched versions that the Sage project generates based
on top of the upstream source.

The only exception to this rule is for *removals* of unused
files or directories.  Some packages contain parts which are not needed
for Sage.  To save space, these may be removed directly from ``src/``.
But be sure to document this in the "Special Update/Build Instructions"
section in ``SPKG.txt``!


.. _section-old-spkg-install:

The File spkg-install
---------------------

The script ``spkg-install`` is run during installation of the Sage
package. In this script, you may make the following assumptions:

- The PATH has the locations of ``sage`` and ``python`` (from the Sage
  installation) at the front. Thus the command::

      python setup.py install

  will run the correct version of Python with everything set up
  correctly. Also, running ``gap`` or ``Singular``, for example, will
  run the correct version.

- The environment variable ``SAGE_ROOT`` points to the root directory
  of the Sage installation.

- The environment variable ``SAGE_LOCAL`` points to the
  ``SAGE_ROOT/local`` directory of the Sage installation.

- The environment variables ``LD_LIBRARY_PATH`` and
  ``DYLD_LIBRARY_PATH`` both have ``SAGE_ROOT/local/lib`` at the
  front.

The ``spkg-install`` script should copy your files to the appropriate
place after doing any build that is necessary.  Here is a template::

    #!/usr/bin/env bash

    if [ -z "$SAGE_LOCAL" ]; then
        echo >&2 "SAGE_LOCAL undefined ... exiting"
        echo >&2 "Maybe run 'sage --sh'?"
        exit 1
    fi

    cd src

    # Apply patches.  See SPKG.txt for information about what each patch
    # does.
    for patch in ../patches/*.patch; do
        [ -r "$patch" ] || continue  # Skip non-existing or non-readable patches
        patch -p1 <"$patch"
        if [ $? -ne 0 ]; then
            echo >&2 "Error applying '$patch'"
            exit 1
        fi
    done

    ./configure --prefix="$SAGE_LOCAL"
    if [ $? -ne 0 ]; then
        echo >&2 "Error configuring PACKAGE_NAME."
        exit 1
    fi

    $MAKE
    if [ $? -ne 0 ]; then
        echo >&2 "Error building PACKAGE_NAME."
        exit 1
    fi

    $MAKE install
    if [ $? -ne 0 ]; then
        echo >&2 "Error installing PACKAGE_NAME."
        exit 1
    fi

    if [ "$SAGE_SPKG_INSTALL_DOCS" = yes ] ; then
        # Before trying to build the documentation, check if any
        # needed programs are present. In the example below, we
        # check for 'latex', but this will depend on the package.
        # Some packages may need no extra tools installed, others
        # may require some.  We use 'command -v' for testing this,
        # and not 'which' since 'which' is not portable, whereas
        # 'command -v' is defined by POSIX.

        # if [ `command -v latex` ] ; then
        #    echo "Good, latex was found, so building the documentation"
        # else
        #    echo "Sorry, can't build the documentation for PACKAGE_NAME as latex is not installed"
        #    exit 1
        # fi


        # make the documentation in a package-specific way
        # for example, we might have
        # cd doc
        # $MAKE html

        if [ $? -ne 0 ]; then
            echo >&2 "Error building PACKAGE_NAME docs."
            exit 1
        fi
        mkdir -p "$SAGE_ROOT/local/share/doc/PACKAGE_NAME"
        # assuming the docs are in doc/*
        cp -R doc/* "$SAGE_ROOT/local/share/doc/PACKAGE_NAME"
    fi


Note that the first line is ``#!/usr/bin/env bash``; this is important
for portability.  Next, the script checks that ``SAGE_LOCAL`` is
defined to make sure that the Sage environment has been set.  After
this, the script may simply run ``cd src`` and then call either
``python setup.py install`` or the autotools sequence
``./configure && make && make install``, or something else along these
lines.

Sometimes, though, it can be more complicated. For example, you might need
to apply the patches from the ``patches`` directory in a particular order. Also,
you should first build (e.g. with ``python setup.py build``, exiting
if there is an error), before installing (e.g. with ``python setup.py
install``). In this way, you would not overwrite a working older
version with a non-working newer version of the spkg.

When copying documentation to
``$SAGE_ROOT/local/share/doc/PACKAGE_NAME``, it may be necessary to
check that only the actual documentation files intended for the user
are copied.  For example, if the documentation is built from ``.tex``
files, you may just need to copy the resulting pdf files, rather than
copying the entire doc directory.  When generating documentation using
Sphinx, copying the ``build/html`` directory generally will copy just
the actual output intended for the user.


.. _section-old-spkg-SPKG-txt:

The File SPKG.txt
-----------------

The old-style ``SPKG.txt`` file is the same as described in
:ref:`section-spkg-SPKG-txt`, but with a hand-maintained changelog
appended since the contents are not part of the Sage repository
tree. It should follow the following pattern::

     == Changelog ==

     Provide a changelog of the spkg here, where the entries have this format:

     === mypackage-0.1.p0 (Mary Smith, 1 Jan 2012) ===

      * Patch src/configure so it builds on Solaris. See Sage trac #137.

     === mypackage-0.1 (Leonhard Euler, 17 September 1783) ===

      * Initial release.  See Sage trac #007.

When the directory (say, ``mypackage-0.1``) is ready, the command

::

    sage --pkg mypackage-0.1

will create the file ``mypackage-0.1.spkg``.  As noted above, this
creates a compressed tar file. Running ``sage --pkg_nc mypackage-0.1``
creates an uncompressed tar file.

When your spkg is ready, you should post about it on ``sage-devel``.
If people there think it is a good idea, then post a link to the spkg
on the Sage trac server (see :ref:`chapter-sage-trac`) so it can be
refereed.  Do not post the spkg itself to the trac server: you only
need to provide a link to your spkg.  If your spkg gets a positive
review, it might be included into the core Sage library, or it might
become an optional download from the Sage website, so anybody can
automatically install it by typing ``sage -i mypackage-version.spkg``.

.. note::

   For any spkg:

   - Make sure that the hg repository contains every file outside the
     ``src`` directory, and that these are all up-to-date and committed
     into the repository.

   - Include an ``spkg-check`` file if possible (see `trac ticket #299`_).

   .. _trac ticket #299: http://trac.sagemath.org/sage_trac/ticket/299

.. note::

    External Magma code goes in ``SAGE_ROOT/src/ext/magma/user``, so
    if you want to redistribute Magma code with Sage as a package that
    Magma-enabled users can use, that is where you would put it. You
    would also want to have relevant Python code to make the Magma
    code easily usable.


.. _section-old-spkg-avoiding-troubles:

Avoiding Troubles
=================

This section contains some guidelines on what an spkg must never do to
a Sage installation. You are encouraged to produce an spkg that is as
self-contained as possible.

#. An spkg must not modify an existing source file in the Sage
   library.
#. Do not allow an spkg to modify another spkg. One spkg can depend on
   other spkg -- see above. You need to first test for the existence of the
   prerequisite spkg before installing an spkg that depends on it.




.. _section-old-spkg-patching-overview:

Overview of Patching SPKGs
==========================

Make sure you are familiar with the structure and conventions relating
to spkg's; see the chapter :ref:`chapter-old-spkg` for
details. Patching an spkg involves patching the installation script of
the spkg and/or patching the upstream source code contained in the
spkg. Say you want to patch the Matplotlib package
``matplotlib-1.0.1.p0``. Note that "p0" denotes the patch level of the
spkg, while "1.0.1" refers to the upstream version of Matplotlib as
contained under ``matplotlib-1.0.1.p0/src/``. The installation script
of that spkg is::

    matplotlib-1.0.1.p0/spkg-install

In general, a script with the name ``spkg-install``  is an
installation script for an spkg. To patch the installation script, use
a text editor to edit that script. Then in the log file ``SPKG.txt``,
provide a high-level description of your changes. Once you are
satisfied with your changes in the installation script and the log
file ``SPKG.txt``, use Mercurial to check in your changes and make
sure to provide a meaningful commit message.

The directory ``src/`` contains the source code provided by the
upstream project. For example, the source code of Matplotlib 1.0.1 is
contained under ::

    matplotlib-1.0.1.p0/src/

To patch the upstream source code, you should edit a copy of the
relevant file -- files in the ``src/`` directory should be untouched,
"vanilla" versions of the source code.  For example, you might copy
the entire ``src/`` directory::

    $ pwd
    matplotlib-1.0.1.p0
    $ cp -pR src src-patched

Then edit files in ``src-patched/``.  Once you are satisfied with your
changes, generate a unified diff between the original file and the
edited one, and save it in ``patches/``::

    $ diff -u src/configure src-patched/configure > patches/configure.patch

Save the unified diff to a file with the same name as the source file
you patched, but using the file extension ".patch". Note that the
directory ``src/`` should not be under revision control, whereas
``patches/`` must be under revision control. The Mercurial
configuration file ``.hgignore`` should contain the following line::

    src/

Ensure that the installation script ``spkg-install`` contains code to
apply the patches to the relevant files under ``src/``. For example,
the file ::

    matplotlib-1.0.1.p0/patches/finance.py.patch

is a patch for the file ::

    matplotlib-1.0.1.p0/src/lib/matplotlib/finance.py

The installation script ``matplotlib-1.0.1.p0/spkg-install`` contains the
following code to install the relevant patches::

    cd src

    # Apply patches.  See SPKG.txt for information about what each patch
    # does.
    for patch in ../patches/*.patch; do
        patch -p1 <"$patch"
        if [ $? -ne 0 ]; then
            echo >&2 "Error applying '$patch'"
            exit 1
        fi
    done

Of course, this could be modified if the order in which the patches
are applied is important, or if some patches were platform-dependent.
For example::

    if [ "$UNAME" = "Darwin" ]; then
        for patch in ../patches/darwin/*.patch; do
            patch -p1 <"$patch"
            if [ $? -ne 0 ]; then
                echo >&2 "Error applying '$patch'"
                exit 1
            fi
        done
    fi

(The environment variable :envvar:`UNAME` is defined by the script
``sage-env``, and is available when ``spkg-install`` is run.)

Now provide a high-level explanation of your changes in ``SPKG.txt``.
Note the format of ``SPKG.txt`` -- see the chapter
:ref:`chapter-old-spkg` for details.  Once you are satisfied with your
changes, use Mercurial to check in your changes with a meaningful
commit message.  Then use the command ``hg tag`` to tag the tip with
the new version number (using "p1" instead of "p0": we have made
changes, so we need to update the patch level)::

    $ hg tag matplotlib-1.0.1.p1

Next, rename the directory ``matplotlib-1.0.1.p0`` to
``matplotlib-1.0.1.p1`` to match the new patch level.  To produce the
actual spkg file, change to the parent directory of
``matplotlib-1.0.1.p1`` and execute ::

    $ /path/to/sage-x.y.z/sage --pkg matplotlib-1.0.1.p1
    Creating Sage package matplotlib-1.0.1.p1

    Created package matplotlib-1.0.1.p1.spkg.

        NAME: matplotlib
     VERSION: 1.0.1.p1
        SIZE: 11.8M
     HG REPO: Good
    SPKG.txt: Good

Spkg files are either bzipped tar files or just plain tar files; the
command ``sage --pkg ...`` produces the bzipped version.  If your spkg
contains mostly binary files which will not compress well, you can use
``sage --pkg_nc ...`` to produce an uncompressed version, i.e., a
plain tar file::

    $ sage --pkg_nc matplotlib-1.0.1.p0/
    Creating Sage package matplotlib-1.0.1.p0/ with no compression

    Created package matplotlib-1.0.1.p0.spkg.

        NAME: matplotlib
     VERSION: 1.0.1.p0
        SIZE: 32.8M
     HG REPO: Good
    SPKG.txt: Good

Note that this is almost three times the size of the compressed
version, so we should use the compressed version!

At this point, you might want to submit your patched spkg for review.
So provide a URL to your spkg on the relevant trac ticket and/or in an
email to the relevant mailing list. Usually, you should not upload
your spkg itself to the relevant trac ticket -- don't post large
binary files to the trac server.


SPKG Versioning
===============

If you want to bump up the version of an spkg, you need to follow some
naming conventions. Use the name and version number as given by the
upstream project, e.g. ``matplotlib-1.0.1``. If the upstream package is
taken from some revision other than a stable version, you need to
append the date at which the revision is made, e.g. the Singular
package ``singular-3-1-0-4-20090818.p3.spkg`` is made with the
revision as of 2009-08-18. If you start afresh from an upstream
release without any patches to its source code, the resulting spkg
need not have any patch-level labels (appending ".p0" is allowed, but
is optional). For example, ``sagenb-0.6.spkg``
is taken from the upstream stable version ``sagenb-0.6`` without any
patches applied to its source code. So you do not see any patch-level
numbering such as ``.p0`` or ``.p1``.

Say you start with ``matplotlib-1.0.1.p0`` and you want to replace
Matplotlib 1.0.1 with version 1.0.2. This entails replacing the source
code for Matplotlib 1.0.1 under ``matplotlib-1.0.1.p0/src/`` with the
new source code. To start with, follow the naming conventions as
described in the section :ref:`section-old-spkg-patching-overview`. If
necessary, remove any obsolete patches and create any new ones,
placing them in the ``patches/`` directory.  Modify the script
``spkg-install`` to take any changes to the patches into account; you
might also have to deal with changes to how the new version of the
source code builds. Then package your replacement spkg using the Sage
command line options ``--pkg`` or ``--pkg_nc`` (or tar and bzip2).

To install your replacement spkg, you use::

    sage -f http://URL/to/package-x.y.z.spkg

or::

    sage -f /path/to/package-x.y.z.spkg

To compile Sage from source with the replacement (standard) spkg,
untar a Sage source tarball, remove the existing spkg under
``SAGE_ROOT/spkg/standard/``. In its place, put your replacement
spkg. Then execute ``make`` from ``SAGE_ROOT``.

