.. comment:
    ***************************************************************************
    If you alter this document, please change the last line:
    **This page was last updated in MONTH YEAR (Sage X.Y).**
    ***************************************************************************

.. HIGHLIGHT:: shell-session

.. _sec-installation-from-sources:

Install from Source Code
========================

.. contents:: Table of contents
   :depth: 2
   :class: this-will-duplicate-information-and-it-is-still-useful-here

Some familiarity with the use of the Unix command line may be required to
build Sage from the :wikipedia:`source code <Source_code>`.

Building Sage from the source code has the major advantage that your install
will be optimized for your particular computer and should therefore offer
better performance and compatibility than a binary install.

Moreover, it offers you full development capabilities:
you can change absolutely any part of Sage or the programs on which it depends,
and recompile the modified parts.

See the file `README.md <https://github.com/sagemath/sage/#readme>`_
in ``SAGE_ROOT`` for information on supported platforms and
step-by-step instructions.

The following sections provide some additional details. Most users
will not need to read them.

.. _section-prereqs:

Prerequisites
-------------

Disk space and memory
^^^^^^^^^^^^^^^^^^^^^

Your computer comes with at least 6 GB of free disk space.
It is recommended to have at least 2 GB of RAM, but you might get away
with less (be sure to have some swap space in this case).

Software prerequisites and recommended packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sage depends on `a large number of software packages
<../reference/spkg/index.html>`_.  Sage provides its own software
distribution providing most of these packages, so you do not have to
worry about having to download and install these packages yourself.

If you extracted Sage from a source tarball, the subdirectory
:file:`upstream` contains the source distributions for all standard
packages on which Sage depends.  If cloned from a git repository, the
upstream tarballs will be downloaded, verified, and cached as part of
the Sage installation process.

However, there are minimal prerequisites for building Sage that
already must be installed on your system:

- `Fundamental system packages required for installing from source
  <../reference/spkg/_prereq.html>`_

- `C/C++ compilers <../reference/spkg/gcc.html>`_

If you have sufficient privileges (for example, on Linux you can
use ``sudo`` to become the ``root`` user), then you can install these packages
using the commands for your platform indicated in the pages linked above.
If you do not have the privileges to do this, ask your system administrator to
do this for you.

In addition to these minimal prerequisites, we strongly recommend to use system
installations of the following:

- `Fortran compiler <../reference/spkg/gfortran.html>`_

- `Python <../reference/spkg/python3.html>`_

Sage developers will also need the `system packages required for
bootstrapping <../reference/spkg/_bootstrap.html>`_; they cannot be
installed by Sage.

When the ``./configure`` script runs, it will check for the presence of many
packages (including the above) and inform you of any that are
missing or have unsuitable versions. **Please read the messages that
``./configure`` prints:** It will inform you which additional system packages
you can install to avoid having to build them from source. This can save a lot of
time.

The following sections provide the commands to install a large
recommended set of packages on various systems, which will minimize
the time it takes to build Sage. This is intended as a convenient
shortcut, but of course you can choose to take a more fine-grained
approach.


.. _sec-installation-from-sources-linux-recommended-installation:

Debian/Ubuntu package installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On Debian ("buster" or newer) or Ubuntu ("bionic" or newer), we recommend that you
install:

.. literalinclude:: debian.txt

If you wish to do Sage development, we recommend that you additionally
install the following:

.. literalinclude:: debian-develop.txt

For all users, we recommend that you install the following system packages,
which provide additional functionality and cannot be installed by Sage:

.. literalinclude:: debian-recommended.txt

In addition to these, if you don't want Sage to build optional packages that might
be available from your OS, cf. the growing list of such packages on :trac:`27330`,
install:

.. literalinclude:: debian-optional.txt

Fedora/Redhat/CentOS package installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On Fedora/Redhat/CentOS, we recommend that you install:

.. literalinclude:: fedora.txt

If you wish to do Sage development, we recommend that you additionally
install the following:

.. literalinclude:: fedora-develop.txt

For all users, we recommend that you install the following system packages,
which provide additional functionality and cannot be installed by Sage:

.. literalinclude:: fedora-recommended.txt

In addition to these, if you don't want Sage to build optional packages that might
be available from your OS, cf. the growing list of such packages on :trac:`27330`,
install:

.. literalinclude:: fedora-optional.txt

Arch Linux package installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On ArchLinux, we recommend that you install:

.. literalinclude:: arch.txt

If you wish to do Sage development, we recommend that you additionally
install the following:

.. literalinclude:: arch-develop.txt

For all users, we recommend that you install the following system packages,
which provide additional functionality and cannot be installed by Sage:

.. literalinclude:: arch-recommended.txt

In addition to these, if you don't want Sage to build optional packages that might
be available from your OS, cf. the growing list of such packages on :trac:`27330`,
install:

.. literalinclude:: arch-optional.txt

OpenSUSE package installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On OpenSUSE, we recommend that you install:

.. literalinclude:: opensuse.txt

If you wish to do Sage development, we recommend that you additionally
install the following:

.. literalinclude:: opensuse-develop.txt

For all users, we recommend that you install the following system packages,
which provide additional functionality and cannot be installed by Sage:

.. literalinclude:: opensuse-recommended.txt

In addition to these, if you don't want Sage to build optional packages that might
be available from your OS, cf. the growing list of such packages on :trac:`27330`,
install:

.. literalinclude:: opensuse-optional.txt

.. _section_macprereqs:

macOS prerequisites
^^^^^^^^^^^^^^^^^^^

On macOS systems, you need a recent version of
`Command Line Tools <https://developer.apple.com/downloads/index.action?=command%20line%20tools>`_.
It provides all the above requirements.

Run the command ``xcode-select --install`` from a Terminal window and click "Install"
in the pop-up dialog box.

If you have already installed `Xcode <https://developer.apple.com/xcode/>`_
(which at the time of writing is freely available in the Mac App Store,
or through https://developer.apple.com/downloads/ provided you registered for an
Apple Developer account), you can install the command line tools from
there as well.

If you have not installed `Xcode <https://developer.apple.com/xcode/>`_
you can get these tools as a relatively small download, but it does require
a registration.

- First, you will need to register as an Apple Developer at
  https://developer.apple.com/register/.

- Having done so, you should be able to download it for free at
  https://developer.apple.com/downloads/index.action?=command%20line%20tools

- Alternately, https://developer.apple.com/opensource/ should have a link
  to Command Line Tools.


macOS package installation
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you use the `Homebrew package manager
<https://brew.sh>`_, you can install the following:

.. literalinclude:: homebrew.txt

Some Homebrew packages are installed "keg-only," meaning that they are
not available in standard paths. To make them accessible when building
Sage, run ::

    $ source SAGE_ROOT/.homebrew-build-env

(replacing ``SAGE_ROOT`` by Sage's home directory). You can add a
command like this to your shell profile if you want the settings to
persist between shell sessions.

If you wish to do Sage development, we recommend that you additionally
install the following:

.. literalinclude:: homebrew-develop.txt

For all users, we recommend that you install the following system packages,
which provide additional functionality and cannot be installed by Sage:

.. literalinclude:: homebrew-recommended.txt

Some additional optional packages are taken care of by:

.. literalinclude:: homebrew-optional.txt


Ubuntu on Windows Subsystem for Linux (WSL) prerequisite installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sage can be installed onto Linux running on Windows Subsystem for Linux (WSL). These instructions describe a fresh install of Ubuntu 20.10, but other distributions or installation methods should work too, though have not been tested.

- Enable hardware-assisted virtualization in the EFI or BIOS of your system. Refer to your system (or motherboard) maker's documentation for instructions on how to do this.

- Set up WSL by following the `official WSL setup guide <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_. Be sure to do the steps to install WSL2 and set it as default.

- Go to the Microsoft Store and install Ubuntu.

- Start Ubuntu from the start menu. Update all packages to the latest version.

- Reboot the all running WSL instances one of the following ways:

  - Open Windows Services and restart the LxssManager service.
  - Open the Command Prompt or Powershell and enter this command::

      wsl --shutdown

- `Upgrade to the Ubuntu 20.10 <https://linuxconfig.org/how-to-upgrade-ubuntu-to-20-10>`_. This step will not be necessary once Ubuntu 20.10 is available in the Microsoft Store.

From this point on, follow the instructions in the :ref:`sec-installation-from-sources-linux-recommended-installation` section.
It is strongly recommended to put the Sage source files in the Linux file system, for example, in the ``/home/username/sage`` directory, and not in the Windows file system (e.g. ``/mnt/c/...``).

You may encounter permission errors of the kind ``"[Errno 13] Permission denied: 'build/bdist.linux-x86_64/wheel/<package>.dist-info'"`` during ``make``.
This usually comes from a permission conflict between the Windows and Linux file system.
To fix it create a temporary build folder in the Linux file system using ``mkdir -p ~/tmp/sage`` and use it for building by ``eval SAGE_BUILD_DIR="~/tmp/sage" make``.
Also see the `related Github issue <https://github.com/pypa/packaging-problems/issues/258>`_ for other workarounds.

When the installation is complete, you may be interested in :ref:`sec-launching-wsl-post-installation`.

.. _section_cygwinprereqs:

Cygwin prerequisite installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Previous versions of Sage targeted the Windows platform using `Cygwin
<https://cygwin.com/>`_.

As of Sage 9.7, we no longer recommend attempting to build Sage on
Cygwin and instead suggest that users on Windows 10 and 11 switch to
installing Sage using Windows Subsystem for Linux (WSL), which gives a
better performance and user/developer experience than Cygwin.

Users on hardware configurations that do not support running WSL, as
well as users on legacy versions of Windows such as Windows 8 may find
it necessary to build Sage on Cygwin.

.. WARNING::

   As of Sage 9.7, :trac:`known issues with several packages
   <query?status=closed&status=needs_info&status=needs_review&status=needs_work&status=new&status=positive_review&component=porting%3A+Cygwin&milestone=sage-9.8&milestone=sage-9.7&milestone=sage-9.6&milestone=sage-9.5&milestone=sage-9.4&milestone=sage-9.3&milestone=sage-9.2&milestone=sage-9.1&col=id&col=summary&col=milestone&col=status&col=priority&col=changetime&col=author&col=reviewer&desc=1&order=changetime>`
   will prevent a successful installation. Users need to be prepared
   to contribute to Sage by fixing these issues.

Use the following instructions to get started.

1.  Download `the 64-bit version of Cygwin <https://cygwin.com/install.html>`_
    (do not get the 32-bit version; it is not supported by Sage).

2.  Run the ``setup-x86_64.exe`` graphical installer.  Pick the default
    options in most cases.  At the package selection screen, use the
    search bar to find and select at least the following packages:
    ``bzip2``, ``coreutils``, ``curl``, ``gawk``, ``gzip``, ``tar``, ``wget``, ``git``.

3.  Start the Cygwin terminal and ensure you get a working bash prompt.

4.  Make sure the path of your Cygwin home directory does not contain
    space characters. Also avoid building in home directories of Windows domain
    users or in paths with capital letters.

    By default, your username in Cygwin is the same as your username in
    Windows.  This might contain spaces and other traditionally
    non-UNIX-friendly characters, e.g., if it is your full name.  You
    can check this as follows::

        $ whoami
        Erik M. Bray

    This means your default home directory on Cygwin contains this
    username verbatim; in the above example, ``/home/Erik M. Bray``.
    It will save some potential trouble if you change your Cygwin home
    directory to contain only alphanumeric characters, for example,
    ``/home/embray``.  The easiest way to do this is to first create
    the home directory you want to use instead, then create an
    ``/etc/passwd`` file specifying that directory as your home, as follows::

        $ whocanibe=embray
        $ mkdir /home/$whocanibe
        $ mkpasswd.exe -l -u "$(whoami)" | sed -r 's,/home/[^:]+,/home/'$whocanibe, > /etc/passwd

    After this, close all Cygwin terminals (ensure nothing in
    ``C:\cygwin64`` is running), then start a new Cygwin terminal and
    your home directory should have moved.

    There are `other ways to do
    this <https://stackoverflow.com/questions/1494658/how-can-i-change-my-cygwin-home-folder-after-installation>`_,
    but the above seems to be the simplest that's still supported.

5.  (Optional) Although it is possible to install Sage's dependencies using the
    Cygwin graphical installer, it is recommended to install the
    `apt-cyg <https://github.com/transcode-open/apt-cyg>`_
    command-line package installer, which is used for the remainder of
    these instructions.  To install ``apt-cyg``, run::

        $ curl -OL https://rawgit.com/transcode-open/apt-cyg/master/apt-cyg
        $ install apt-cyg /usr/local/bin
        $ rm -f apt-cyg

6.  Then, to install the current set of system packages known to work for building
    Sage, run the following command (or use the graphical installer to
    select and install these packages):

    .. literalinclude:: cygwin.txt

    Optional packages that are also known to be installable via system packages
    include:

    .. literalinclude:: cygwin-optional.txt

.. NOTE::

   On Cygwin, at any point in time after building/installing software,
   it may be required to  "rebase" ``dll`` files.
   Sage provides some scripts, located in :file:`$SAGE_LOCAL/bin`, to do so:

   - ``sage-rebaseall.sh``, a shell script which calls Cygwin's
     ``rebaseall`` program.  It must be run within a ``dash`` shell
     from the :envvar:`SAGE_ROOT` directory after all other Cygwin
     processes have been shut down and needs write-access to the
     system-wide rebase database located at
     :file:`/etc/rebase.db.i386`, which usually means administrator
     privileges.  It updates the system-wide database and adds Sage
     dlls to it, so that subsequent calls to ``rebaseall`` will take
     them into account.

   - ``sage-rebase.sh``, a shell script which calls Cygwin's ``rebase`` program
     together with the ``-O/--oblivious`` option.
     It must be run within a shell from :envvar:`SAGE_ROOT` directory.
     Contrary to the ``sage-rebaseall.sh`` script, it neither updates the
     system-wide database, nor adds Sage dlls to it.
     Therefore, subsequent calls to ``rebaseall`` will not take them into account.

   - ``sage-rebaseall.bat`` (respectively ``sage-rebase.bat``), an MS-DOS batch
     file which calls the ``sage-rebaseall.sh`` (respectively ``sage-rebase.sh``)
     script.
     It must be run from a Windows command prompt, after adjusting
     :envvar:`SAGE_ROOT` to the Windows location of Sage's home directory, and, if
     Cygwin is installed in a non-standard location, adjusting
     :envvar:`CYGWIN_ROOT` as well.

   Some systems may encounter this problem frequently enough to make building or
   testing difficult.
   If executing the above scripts or directly calling ``rebaseall`` does not solve
   rebasing issues, deleting the system-wide database and then regenerating it
   from scratch, e.g., by executing ``sage-rebaseall.sh``, might help.


Other platforms
^^^^^^^^^^^^^^^

On Solaris, you would use ``pkgadd`` and on OpenSolaris ``ipf`` to install
the necessary software.

On other systems, check the documentation for your particular operating system.

.. _section_conda_compilers:



Notes on using conda
^^^^^^^^^^^^^^^^^^^^

If you don't want conda to be used by sage, deactivate conda (for the current shell session).

  - Type::

      $ conda deactivate

  - Repeat the command until ``conda info`` shows::

      $ conda info

      active environment : None
      ...

  Then SageMath will be built either using the compilers provided by the
  operating system, or its own compilers.


Additional software
-------------------

Recommended programs
^^^^^^^^^^^^^^^^^^^^

The following programs are recommended.
They are not strictly required at build time or at run time,
but provide additional capabilities:

- **dvipng**.
- **ffmpeg**.
- **ImageMagick**.
- **LaTeX**: highly recommended.

It is highly recommended that you have
:wikipedia:`LaTeX <LaTeX>`
installed, but it is not required.
The most popular packaging is `TeX Live <https://www.tug.org/texlive/>`_,
which can be installed following the directions on their web site.
On Linux systems you can alternatively install your distribution's
texlive packages::

    $ sudo apt-get install texlive       # debian
    $ sudo yum install texlive           # redhat

or similar commands. In addition to the base TeX Live install, you may
need some optional TeX Live packages, for example
country-specific babel packages for the localized Sage
documentation.

If you don't have either ImageMagick or ffmpeg, you won't be able to
view animations.
ffmpeg can produce animations in more different formats than ImageMagick,
and seems to be faster than ImageMagick when creating animated GIFs.
Either ImageMagick or dvipng is used for displaying some LaTeX output in the
Sage notebook.

On Debian/Ubuntu, the following system packages are recommended.

- ``texlive-generic-extra`` (to generate pdf documentation)

- ``texlive-xetex`` (to convert Jupyter notebooks to pdf)

- ``latexmk`` (to generate pdf documentation)

- ``pandoc`` (to convert Jupyter notebooks to pdf)

- ``dvipng`` (to render text with LaTeX in Matplotlib)

- ``default-jdk`` (to run the Jmol 3D viewer from the console and generate images for 3D plots in the documentation)

- ``ffmpeg`` (to produce animations)

- ``libavdevice-dev`` (to produce animations)

Tcl/Tk
^^^^^^

If you want to use `Tcl/Tk <https://www.tcl.tk/>`_ libraries in Sage,
you need to install the Tcl/Tk and its development headers before building
Sage.
Sage's Python will then automatically recognize your system's install of
Tcl/Tk.

On Linux systems, these are usually provided by the **tk** and **tk-dev**
(or **tk-devel**) packages which can be installed using::

    $ sudo apt-get install tk tk-dev

or similar commands.

If you installed Sage first, all is not lost. You just need to rebuild
Sage's Python and any part of Sage relying on it::

    $ sage -f python3  # rebuild Python3
    $ make             # rebuild components of Sage depending on Python

after installing the Tcl/Tk development libraries as above.

If

.. skip

.. CODE-BLOCK:: ipycon

   sage: import _tkinter
   sage: import Tkinter

does not raise an ``ImportError``, then it worked.

.. _build-from-source-step-by-step:

Step-by-step installation procedure
-----------------------------------

General procedure
^^^^^^^^^^^^^^^^^

#. Follow the procedure in the file `README.md <https://github.com/sagemath/sage/#readme>`_
   in ``SAGE_ROOT``.

#. Additional remarks:
   You do not need to be logged in as root, since no files are
   changed outside of the :file:`sage-x.y` directory.
   In fact, **it is inadvisable to build Sage as root**, as the root account
   should only be used when absolutely necessary and mistyped commands can have
   serious consequences if you are logged in as root.

   Typing ``make`` performs the usual steps for each Sage's dependency,
   but installs all the resulting files into the installation prefix.
   Depending on the age and the architecture of your system, it can take from
   a few tens of minutes to several hours to build Sage from source.
   On really slow hardware, it can even take a few days to build Sage.

   Each component of Sage has its own build log, saved in
   :file:`SAGE_ROOT/logs/pkgs`.
   If the build of Sage fails, you will see a message mentioning which
   package(s) failed to build and the location of the log file for each
   failed package.
   If this happens, then paste the contents of these log file(s)
   to the Sage support
   newsgroup at https://groups.google.com/group/sage-support.
   If the log files are very large (and many are), then don't paste the whole
   file, but make sure to include any error messages.
   It would also be helpful to include the type of operating system
   (Linux, macOS, Solaris, OpenSolaris, Cygwin, or any other system),
   the version and release date of that operating system and the version of
   the copy of Sage you are using.
   (There are no formal requirements for bug reports -- just send them;
   we appreciate everything.)

   See :ref:`section_make` for some targets for the ``make`` command and
   :ref:`section_envvar` for additional information on useful environment
   variables used by Sage.

#. To start Sage, you can now simply type from Sage's home directory::

       $ ./sage

   You should see the Sage prompt, which will look something like this::

       $ sage
       ┌────────────────────────────────────────────────────────────────────┐
       │ SageMath version 8.8, Release Date: 2019-06-26                     │
       │ Using Python 3.10.4. Type "help()" for help.                       │
       └────────────────────────────────────────────────────────────────────┘
       sage:

   Note that Sage should take well under a minute when it starts for the first
   time, but can take several minutes if the file system is slow or busy.
   Since Sage opens a lot of files, it is preferable to install Sage on a fast
   filesystem if possible.

   Just starting successfully tests that many of the components built
   correctly.
   Note that this should have been already automatically tested during the
   build process.
   If the above is not displayed (e.g., if you get a massive traceback), please
   report the problem, e.g., at https://groups.google.com/group/sage-support.

   After Sage has started, try a simple command:

   .. CODE-BLOCK:: ipycon

       sage: 2 + 2
       4

   Or something slightly more complicated:

   .. CODE-BLOCK:: ipycon

       sage: factor(2005)
       5 * 401


#. Optional, but highly recommended:
   Test the install by typing ``./sage --testall``.
   This runs most examples in the source code and makes sure that they run
   exactly as claimed.
   To test all examples, use ``./sage --testall --optional=all --long``;
   this will run examples that take a long time, and those that depend on
   optional packages and software, e.g., Mathematica or Magma.
   Some (optional) examples will therefore likely fail.

   Alternatively, from within :file:`$SAGE_ROOT`, you can type ``make test``
   (respectively ``make ptest``) to run all the standard test code serially
   (respectively in parallel).

   Testing the Sage library can take from half an hour to several hours,
   depending on your hardware.
   On slow hardware building and testing Sage can even take several days!


#. Optional:
   Check the interfaces to any other software that you have available.
   Note that each interface calls its corresponding program by a particular
   name: `Mathematica <https://www.wolfram.com/mathematica/>`_ is invoked by
   calling ``math``, `Maple <https://www.maplesoft.com/>`_ by calling ``maple``,
   etc.
   The easiest way to change this name or perform other customizations is
   to create a redirection script in :file:`$SAGE_ROOT/local/bin`.
   Sage inserts this directory at the front of your :envvar:`PATH`, so your
   script may need to use an absolute path to avoid calling itself; also, your
   script should pass along all of its arguments.
   For example, a ``maple`` script might look like:

   .. CODE-BLOCK:: bash

       #!/bin/sh

       exec /etc/maple10.2/maple.tty "$@"

#. Optional:
   There are different possibilities to make using Sage a little easier:

   - Make a symbolic link from :file:`/usr/local/bin/sage` (or another
     directory in your :envvar:`PATH`) to :file:`$SAGE_ROOT/sage`::

         $ ln -s /path/to/sage-x.y/sage /usr/local/bin/sage

     Now simply typing ``sage`` from any directory should be sufficient to run
     Sage.

   - Copy :file:`$SAGE_ROOT/sage` to a location in your :envvar:`PATH`.
     If you do this, make sure you edit the line:

     .. CODE-BLOCK:: bash

         #SAGE_ROOT=/path/to/sage-version

     at the beginning of the copied ``sage`` script according to the direction
     given there to something like:

     .. CODE-BLOCK:: bash

         SAGE_ROOT=<SAGE_ROOT>

     (note that you have to change ``<SAGE_ROOT>`` above!).
     It is best to edit only the copy, not the original.

   - For `KDE <https://www.kde.org/>`_ users, create a bash script called
     ``sage`` containing the lines
     (note that you have to change ``<SAGE_ROOT>`` below!):

     .. CODE-BLOCK:: bash

         #!/usr/bin/env bash

         konsole -T "sage" -e <SAGE_ROOT>/sage

     make it executable::

         $ chmod a+x sage

     and put it somewhere in your :envvar:`PATH`.

     You can also make a KDE desktop icon with this line as the command
     (under the Application tab of the Properties of the icon, which you get my
     right clicking the mouse on the icon).

   - On Linux and macOS systems, you can make an alias to
     :file:`$SAGE_ROOT/sage`.
     For example, put something similar to the following line in your
     :file:`.bashrc` file:

     .. CODE-BLOCK:: bash

         alias sage=<SAGE_ROOT>/sage

     (Note that you have to change ``<SAGE_ROOT>`` above!)
     Having done so, quit your terminal emulator and restart it.
     Now typing ``sage`` within your terminal emulator should start Sage.

#. Optional:
   Install optional Sage packages and databases. See `the list of optional packages
   in the reference manual <../reference/spkg/index.html#optional-packages>`_ for
   detailed information, or type ``sage --optional`` (this requires an Internet connection).

   Then type ``sage -i <package-name>`` to automatically download and install
   a given package.

#. Have fun! Discover some amazing conjectures!


.. _section_make:

Make targets
------------

To build Sage from scratch, you would typically execute ``make`` in Sage's home
directory to build Sage and its :wikipedia:`HTML <HTML>`
documentation.
The ``make`` command is pretty smart, so if your build of Sage is interrupted,
then running ``make`` again should cause it to pick up where it left off.
The ``make`` command can also be given options, which control what is built and
how it is built:

- ``make build`` builds Sage: it compiles all of the Sage packages.
  It does not build the documentation.

- ``make doc`` builds Sage's documentation in HTML format.
  Note that this requires that Sage be built first, so it will automatically
  run ``make build`` first.
  Thus, running ``make doc`` is equivalent to running ``make``.

- ``make doc-pdf`` builds Sage's documentation in PDF format. This also
  requires that Sage be built first, so it will automatically run ``make
  build``.

- ``make doc-html-no-plot`` builds Sage's documentation in html format
  but skips the inclusion of graphics auto-generated using the
  ``.. PLOT`` markup and the ``sphinx_plot`` function. This is
  primarily intended for use when producing certain binary
  distributions of Sage, to lower the size of the distribution. As of
  this writing (December 2014, Sage 6.5), there are only a few such
  plots, adding about 4M to the :file:`local/share/doc/sage/` directory.
  In the future, this may grow, of course. Note: after using this, if you
  want to build the documentation and include the pictures, you should
  run ``make doc-uninstall``, because the presence, or lack, of pictures
  is cached in the documentation output.
  You can benefit from this no-plot feature with other make targets by doing
  ``export SAGE_DOCBUILD_OPTS+=' --no-plot'``

- ``make ptest`` and ``make ptestlong``: these run Sage's test suite.
  The first version skips tests that need more than a few seconds to complete
  and those which depend on optional packages or additional software.
  The second version includes the former, and so it takes longer.
  The "p" in ``ptest`` stands for "parallel": tests are run in parallel.
  If you want to run tests serially, you can use ``make test`` or
  ``make testlong`` instead.
  If you want to run tests depending on optional packages and additional
  software, you can use ``make testall``, ``make ptestall``,
  ``make testalllong``, or ``make ptestalllong``.

- ``make doc-uninstall`` and ``make doc-clean`` each remove several
  directories which are produced when building the documentation.

- ``make distclean`` restores the Sage directory to its state before doing any
  building: it is almost equivalent to deleting Sage's entire home directory and
  unpacking the source tarfile again, the only difference being that the
  :file:`.git` directory is preserved, so git branches are not deleted.

.. _section_envvar:

Environment variables
---------------------

Sage uses several environment variables to control its build process.
Most users won't need to set any of these: the build process just works on many
platforms.
(Note though that setting :envvar:`MAKE`, as described below, can significantly
speed up the process.)
Building Sage involves building about 100 packages, each of which has its own
compilation instructions.

The Sage source tarball already includes the sources for all standard
packages, that is, it allows you to build Sage without internet
connection. The git repository, however, does not contain the source
code for third-party packages. Instead, it will be downloaded as
needed (Note: you can run ``make download`` to force downloading
packages before building). Package downloads use the Sage mirror
network, the nearest mirror will be determined automatically for
you. This is influenced by the following environment variable:

- :envvar:`SAGE_SERVER` - Try the specified mirror first, before
  falling back to the official Sage mirror list. Note that Sage will
  search the directory

  - ``SAGE_SERVER/spkg/upstream``

  for upstream tarballs.

Here are some of the more commonly used variables affecting the build process:

- :envvar:`MAKE` - one useful setting for this variable when building Sage is
  ``MAKE='make -jNUM'`` to tell the ``make`` program to run ``NUM`` jobs in
  parallel when building.
  Note that some Sage packages may not support this variable.

  Some people advise using more jobs than there are CPU cores, at least if the
  system is not heavily loaded and has plenty of RAM; for example, a good
  setting for ``NUM`` might be between 1 and 1.5 times the number of cores.
  In addition, the ``-l`` option sets a load limit: ``MAKE='make -j4 -l5.5``,
  for example, tells ``make`` to try to use four jobs, but to not start more
  than one job if the system load average is above 5.5.
  See the manual page for GNU ``make``: `Command-line options
  <https://www.gnu.org/software/make/manual/make.html#Options-Summary>`_
  and `Parallel building
  <https://www.gnu.org/software/make/manual/make.html#Parallel>`_.

  .. warning::

      Some users on single-core macOS machines have reported problems when
      building Sage with ``MAKE='make -jNUM'`` with ``NUM`` greater than one.

- :envvar:`SAGE_NUM_THREADS` - if set to a number, then when building the
  documentation, parallel doctesting, or running ``sage -b``, use this many
  threads.
  If this is not set, then determine the number of threads using the value of
  the :envvar:`MAKE` (see above) or :envvar:`MAKEFLAGS` environment variables.
  If none of these specifies a number of jobs, use one thread (except for
  parallel testing: there we use a default of the number of CPU cores, with a
  maximum of 8 and a minimum of 2).

- :envvar:`V` - if set to ``0``, silence the build.  Instead of
  showing a detailed compilation log, only one line of output is shown
  at the beginning and at the end of the installation of each Sage
  package.  To see even less output, use::

    $ make -s V=0

  (Note that the above uses the syntax of setting a Makefile variable.)

- :envvar:`SAGE_CHECK` - if set to ``yes``, then during the build process,
  or when installing packages manually,
  run the test suite for each package which has one, and stop with an error
  if tests are failing.  If set to ``warn``, then only a warning is printed
  in this case.
  See also :envvar:`SAGE_CHECK_PACKAGES`.

- :envvar:`SAGE_CHECK_PACKAGES` - if :envvar:`SAGE_CHECK` is set to ``yes``,
  then the default behavior is to run test suites for all spkgs which contain
  them.
  If :envvar:`SAGE_CHECK_PACKAGES` is set, it should be a comma-separated list
  of strings of the form ``package-name`` or ``!package-name``.
  An entry ``package-name`` means to run the test suite for the named package
  regardless of the setting of :envvar:`SAGE_CHECK`.
  An entry ``!package-name`` means to skip its test suite.
  So if this is set to ``ppl,!python3``, then always run the test suite for
  PPL, but always skip the test suite for Python 3.

  .. note::

     As of Sage 9.1, the test suites for the Python 2 and 3 spkgs fail
     on most platforms.  So when this variable is empty or unset, Sage
     uses a default of ``!python2,!python3``.

- :envvar:`SAGE_INSTALL_GCC` - **Obsolete, do not use, to be removed**

- :envvar:`SAGE_INSTALL_CCACHE` - by default Sage doesn't install ccache,
  however by setting ``SAGE_INSTALL_CCACHE=yes`` Sage will install ccache.
  Because the Sage distribution is quite large, the maximum cache is set to 4G.
  This can be changed by running ``sage -sh -c "ccache --max-size=SIZE"``,
  where ``SIZE`` is specified in gigabytes, megabytes, or kilobytes by
  appending a "G", "M", or "K".

  Sage does not include the sources for ccache since it is an optional package.
  Because of this, it is necessary to have an Internet connection while
  building ccache for Sage, so that Sage can pull down the necessary
  sources.

- .. _sage_debug:

  :envvar:`SAGE_DEBUG` - controls debugging support.
  There are three different possible values:

  * Not set (or set to anything else than "yes" or "no"): build binaries with
    debugging symbols, but no special debug builds.
    This is the default.
    There is no performance impact, only additional disk space is used.

  * ``SAGE_DEBUG=no``: ``no`` means no debugging symbols (that is, no
    ``gcc -g``), which saves some disk space.

  * ``SAGE_DEBUG=yes``: build debug versions if possible (in particular,
    Python is built with additional debugging turned on and Singular is built
    with a different memory manager).
    These will be notably slower but, for example, make it much easier to
    pinpoint memory allocation problems.

  Instead of using :envvar:`SAGE_DEBUG` one can configure with
  ``--enable-debug={no|symbols|yes}``.

- :envvar:`SAGE_PROFILE` - controls profiling support. If this is set
  to ``yes``, profiling support is enabled where possible. Note that
  Python-level profiling is always available; This option enables
  profiling in Cython modules.

- :envvar:`SAGE_SPKG_INSTALL_DOCS` - if set to ``yes``, then install
  package-specific documentation to
  :file:`$SAGE_ROOT/local/share/doc/PACKAGE_NAME/` when an spkg is
  installed.
  This option may not be supported by all spkgs.
  Some spkgs might also assume that certain programs are available on the
  system (for example, ``latex`` or ``pdflatex``).

- :envvar:`SAGE_DOCBUILD_OPTS` - the value of this variable is passed as an
  argument to ``sage --docbuild all html`` or ``sage --docbuild all pdf`` when
  you run ``make``, ``make doc``, or ``make doc-pdf``.
  For example, you can add ``--no-plot`` to this variable to avoid building
  the graphics coming from the ``.. PLOT`` directive within the documentation,
  or you can add ``--include-tests-blocks`` to include all "TESTS" blocks in the
  reference manual. Run ``sage --docbuild help`` to see the full list
  of options.

- :envvar:`SAGE_BUILD_DIR` - the default behavior is to build each spkg in a
  subdirectory of :file:`$SAGE_ROOT/local/var/tmp/sage/build/`; for
  example, build version 7.27.0 of
  :file:`ipython` in the directory
  :file:`$SAGE_ROOT/local/var/tmp/sage/build/ipython-7.27.0/`.
  If this variable is set, then build in
  :file:`$SAGE_BUILD_DIR/ipython-7.27.0/` instead.
  If the directory :file:`$SAGE_BUILD_DIR` does not exist, it is created.
  As of this writing (Sage 4.8), when building the standard Sage packages,
  1.5 gigabytes of free space are required in this directory (or more if
  ``SAGE_KEEP_BUILT_SPKGS=yes`` -- see below); the exact amount of required
  space varies from platform to platform.
  For example, the block size of the file system will affect the amount of
  space used, since some spkgs contain many small files.

  .. warning::

      The variable :envvar:`SAGE_BUILD_DIR` must be set to the full path name
      of either an existing directory for which the user has write permissions,
      or to the full path name of a nonexistent directory which the user has
      permission to create.
      The path name must contain **no spaces**.

- :envvar:`SAGE_KEEP_BUILT_SPKGS` - the default behavior is to delete each
  build directory -- the appropriate subdirectory of
  :file:`$SAGE_ROOT/local/var/tmp/sage/build` or
  :file:`$SAGE_BUILD_DIR` -- after each spkg
  is successfully built, and to keep it if there were errors installing the
  spkg.
  Set this variable to ``yes`` to keep the subdirectory regardless.
  Furthermore, if you install an spkg for which there is already a
  corresponding subdirectory, for example left over from a previous build,
  then the default behavior is to delete that old subdirectory.
  If this variable is set to ``yes``, then the old subdirectory is moved to
  :file:`$SAGE_ROOT/local/var/tmp/sage/build/old/` (or
  :file:`$SAGE_BUILD_DIR/old`),
  overwriting any already existing file or directory with the same name.

  .. note::

      After a full build of Sage (as of version 4.8), these subdirectories can
      take up to 6 gigabytes of storage, in total, depending on the platform
      and the block size of the file system.
      If you always set this variable to ``yes``, it can take even more space:
      rebuilding every spkg would use double the amount of space, and any
      upgrades to spkgs would create still more directories, using still more
      space.

  .. note::

      In an existing Sage installation, running ``sage -i -s <package-name>``
      or ``sage -f -s <package-name>`` installs the spkg ``<package-name>`` and
      keeps the corresponding build directory; thus setting
      :envvar:`SAGE_KEEP_BUILT_SPKGS` to ``yes`` mimics this behavior when
      building Sage from scratch or when installing individual spkgs.
      So you can set this variable to ``yes`` instead of using the ``-s`` flag
      for ``sage -i`` and ``sage -f``.

- .. _sage_fat_binary:

  :envvar:`SAGE_FAT_BINARY` - to build binaries that will run on the
  widest range of target CPUs set this variable to ``yes`` before
  building Sage or configure with ``--enable-fat-binary``.
  This does not make the binaries relocatable, it only
  avoids newer CPU instruction set extensions. For relocatable (=can
  be moved to a different directory) binaries, you must use
  https://github.com/sagemath/binary-pkg

- :envvar:`SAGE_SUDO` - set this to ``sudo -E`` or to any other
  command prefix that is necessary to write into a installation
  hierarchy (:envvar:`SAGE_LOCAL`) owned by root or another user.
  Note that this command needs to preserve environment variable
  settings (plain ``sudo`` does not).

  Not all Sage packages currently support :envvar:`SAGE_SUDO`.

  Therefore this environment variable is most useful when a system
  administrator wishes to install an additional Sage package that
  supports :envvar:`SAGE_SUDO`, into a root-owned installation
  hierarchy (:envvar:`SAGE_LOCAL`).

Environment variables dealing with specific Sage packages:

- :envvar:`SAGE_MATPLOTLIB_GUI` - if set to anything non-empty except ``no``,
  then Sage will attempt to build the graphical backend when it builds the
  matplotlib package.

- :envvar:`PARI_CONFIGURE` - use this to pass extra parameters to
  PARI's ``Configure`` script, for example to specify graphics
  support (which is disabled by default). See the file
  :file:`build/pkgs/pari/spkg-install` for more information.

- :envvar:`SAGE_TUNE_PARI`: If yes, enable PARI self-tuning. Note that
  this can be time-consuming. If you set this variable to "yes", you
  will also see this: ``WARNING: Tuning PARI/GP is unreliable. You may
  find your build of PARI fails, or PARI/GP does not work properly
  once built. We recommend to build this package with
  SAGE_CHECK="yes".``

- :envvar:`PARI_MAKEFLAGS`: The value of this variable is passed as an
  argument to the ``$MAKE`` command when compiling PARI.

Some standard environment variables which are used by Sage:

- :envvar:`CC` - while some programs allow you to use this to specify your C
  compiler, **not every Sage package recognizes this**.
  If GCC is installed within Sage, :envvar:`CC` is ignored and Sage's ``gcc``
  is used instead.

- :envvar:`CPP` - similarly, this will set the C preprocessor for some Sage
  packages, and similarly, using it is likely quite risky.
  If GCC is installed within Sage, :envvar:`CPP` is ignored and Sage's ``cpp``
  is used instead.

- :envvar:`CXX` - similarly, this will set the C++ compiler for some Sage
  packages, and similarly, using it is likely quite risky.
  If GCC is installed within Sage, :envvar:`CXX` is ignored and Sage's ``g++``
  is used instead.

- :envvar:`FC` - similarly, this will set the Fortran compiler.
  This is supported by all Sage packages which have Fortran code.
  However, for historical reasons, the value is hardcoded during the initial
  ``make`` and subsequent changes to ``$FC`` might be ignored (in which case,
  the original value will be used instead).
  If GCC is installed within Sage, :envvar:`FC` is ignored and Sage's
  ``gfortran`` is used instead.

- :envvar:`CFLAGS`, :envvar:`CXXFLAGS` and :envvar:`FCFLAGS` - the flags for
  the C compiler, the C++ compiler and the Fortran compiler, respectively.
  The same comments apply to these: setting them may cause problems, because
  they are not universally respected among the Sage packages. Note
  also that ``export CFLAGS=""`` does not have the same effect as
  ``unset CFLAGS``. The latter is preferable.

- Similar comments apply to other compiler and linker flags like
  :envvar:`CPPFLAGS`, :envvar:`LDFLAGS`, :envvar:`CXXFLAG64`,
  :envvar:`LDFLAG64`, and :envvar:`LD`.

- :envvar:`OPENBLAS_CONFIGURE` - adds additional configuration flags for
  the OpenBLAS package that gets added to the make command. (see :trac:`23272`)

Variables dealing with doctesting:

- :envvar:`SAGE_TIMEOUT` - used for Sage's doctesting: the number of seconds
  to allow a doctest before timing it out.
  If this isn't set, the default is 300 seconds (5 minutes).

- :envvar:`SAGE_TIMEOUT_LONG` - used for Sage's doctesting: the number of
  seconds to allow a doctest before timing it out, if tests are run using
  ``sage -t --long``.
  If this isn't set, the default is 1800 seconds (30 minutes).

- :envvar:`SAGE_TEST_GLOBAL_ITER`, :envvar:`SAGE_TEST_ITER`: these can
  be used instead of passing the flags ``--global-iterations`` and
  ``--file-iterations``, respectively, to ``sage -t``. Indeed, these
  variables are only used if the flags are unset. Run ``sage -t -h``
  for more information on the effects of these flags (and therefore
  these variables).

Sage sets some other environment variables. The most accurate way to
see what Sage does is to first run ``env`` from a shell prompt to see
what environment variables you have set. Then run ``sage --sh -c
env`` to see the list after Sage sets its variables. (This runs a
separate shell, executes the shell command ``env``, and then exits
that shell, so after running this, your settings will be restored.)
Alternatively, you can peruse the shell script
:file:`src/bin/sage-env`.

Sage also has some environment-like settings. Some of these correspond
to actual environment variables while others have names like
environment variables but are only available while Sage is running. To
see a list, execute ``sage.env.[TAB]`` while running Sage.

.. comment:
    ***************************************************************************
    FIX THIS!

    Variables dealing with valgrind and friends:

    - :envvar:`SAGE_TIMEOUT_VALGRIND` - used for Sage's doctesting: the
      number of seconds to allow a doctest before timing it out, if tests
      are run using ``??``.  If this isn't set, the default is 1024*1024
      seconds.

    - :envvar:`SAGE_VALGRIND` - trigger black magic in Python.

    - :envvar:`SAGE_MEMCHECK_FLAGS`, :envvar:`SAGE_MASSIF_FLAGS`,
      :envvar:`SAGE_CACHEGRIND_FLAGS`, :envvar:`SAGE_OMEGA_FLAGS` - flags
      used when using valgrind and one of the tools "memcheck", "massif",
      "cachegrind", or "omega"
    ***************************************************************************


Installation in a Multiuser Environment
---------------------------------------

This section addresses the question of how a system administrator can install
a single copy of Sage in a multi-user computer network.

#. Using ``sudo``, create the installation directory, for example,
   ``/opt/sage/sage-x.y``. We refer to it as ``SAGE_LOCAL`` in the
   instructions below. Do not try to install into a directory that
   already contains other software, such as ``/usr/local``::

       $ sudo mkdir -p SAGE_LOCAL

#. Make the directory writable for you and readable by everyone::

       $ sudo chown $(id -un) SAGE_LOCAL
       $ sudo chmod 755 SAGE_LOCAL

#. Build and install Sage, following the instructions in `README.md
   <https://github.com/sagemath/sage/#readme>`_, using the
   ``configure`` option ``--prefix=SAGE_LOCAL``.

   Do not use ``sudo`` for this step; building Sage must be done using
   your normal user account.

#. Optionally, create a symbolic link to the installed ``sage`` script
   in a directory that is in the users' :envvar:`PATH`, for example
   ``/usr/local/bin``::

       $ sudo ln -s SAGE_LOCAL/bin/sage /usr/local/bin/sage

#. Optionally, change permissions to prevent accidental changes to
   the installation by yourself::

       $ sudo chown -R root SAGE_LOCAL


**This page was last updated in September 2022 (Sage 9.8).**
