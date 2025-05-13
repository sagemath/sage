.. HIGHLIGHT:: shell-session

.. _sec-installation-from-sources:

Install from Source Code
========================

Building Sage from the source code has the major
advantage that your install will be optimized for your particular computer and
should therefore offer better performance and compatibility than a binary
install.

Moreover, it offers you full development capabilities: you can change
absolutely any part of Sage or the packages on which it depends, and recompile
the modified parts.

See the file `README.md <https://github.com/sagemath/sage/#readme>`_
in ``SAGE_ROOT`` for information on supported platforms and
step-by-step instructions.

The following sections provide some additional details. Most users will not
need to read them. Some familiarity with the use of the Unix command line may
be required to build Sage from the source code.


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
./configure prints:** It will inform you which additional system packages
you can install to avoid having to build them from source. This can save a lot of
time.

The following sections provide the commands to install a large
recommended set of packages on various systems, which will minimize
the time it takes to build Sage. This is intended as a convenient
shortcut, but of course you can choose to take a more fine-grained
approach.


.. _sec-installation-from-sources-linux-recommended-installation:

Linux system package installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend that you install the following packages, depending on your distribution:

.. tab:: Debian/Ubuntu

   .. literalinclude:: debian.txt

.. tab:: Fedora/Redhat/CentOS

   .. literalinclude:: fedora.txt

.. tab:: Arch Linux

   .. literalinclude:: arch.txt

.. tab:: OpenSUSE

   .. literalinclude:: opensuse.txt

.. tab:: Void Linux

   .. literalinclude:: void.txt

If you wish to do Sage development, we recommend that you additionally
install the following:

.. tab:: Debian/Ubuntu

   .. literalinclude:: debian-develop.txt

.. tab:: Fedora/Redhat/CentOS

   .. literalinclude:: fedora-develop.txt

.. tab:: Arch Linux

   .. literalinclude:: arch-develop.txt

.. tab:: OpenSUSE

   .. literalinclude:: opensuse-develop.txt

.. tab:: Void Linux

   .. literalinclude:: void-develop.txt

For all users, we recommend that you install the following system
packages, which provide additional functionality and cannot be
installed by Sage.  In particular, this includes :wikipedia:`LaTeX
<LaTeX>` and related tools. In addition to a base install of :ref:`TeX
Live <spkg_texlive>`, our lists of system packages below include
everything that is needed for generating the Sage documentation in PDF
format.  For converting Jupyter notebooks to PDF, also the document
converter :ref:`pandoc <spkg_pandoc>` is needed.  For making
animations, Sage needs to use one of the packages :ref:`FFmpeg
<spkg_ffmpeg>` and :ref:`ImageMagick <spkg_imagemagick>`.

.. tab:: Debian/Ubuntu

   .. literalinclude:: debian-recommended.txt

.. tab:: Fedora/Redhat/CentOS

   .. literalinclude:: fedora-recommended.txt

.. tab:: Arch Linux

   .. literalinclude:: arch-recommended.txt

.. tab:: OpenSUSE

   .. literalinclude:: opensuse-recommended.txt

.. tab:: Void Linux

   .. literalinclude:: void-recommended.txt

In addition to these, if you don't want Sage to build optional packages that might
be available from your OS, cf. the growing list of such packages on :issue:`27330`,
install:

.. tab:: Debian/Ubuntu

   .. literalinclude:: debian-optional.txt

.. tab:: Fedora/Redhat/CentOS

   .. literalinclude:: fedora-optional.txt

.. tab:: Arch Linux

   .. literalinclude:: arch-optional.txt

.. tab:: OpenSUSE

   .. literalinclude:: opensuse-optional.txt

.. tab:: Void Linux

   .. literalinclude:: void-optional.txt


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

WSL prerequisites
^^^^^^^^^^^^^^^^^

Ubuntu on Windows Subsystem for Linux (WSL) prerequisite installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Refer to :ref:`installation-guide` for installing Ubuntu on
Windows Subsystem for Linux (WSL). These instructions describe a fresh
install of Ubuntu, the default distribution in WSL, but other
distributions or installation methods should work too.

From this point on, follow the instructions in the :ref:`sec-installation-from-sources-linux-recommended-installation` section.
It is strongly recommended to put the Sage source files in the Linux file system, for example, in the ``/home/username/sage`` directory, and not in the Windows file system (e.g. ``/mnt/c/...``).

WSL permission denied error when building ``packaging`` package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may encounter permission errors of the kind ``"[Errno 13] Permission denied: 'build/bdist.linux-x86_64/wheel/<package>.dist-info'"`` during ``make``.
This usually comes from a permission conflict between the Windows and Linux file system.
To fix it create a temporary build folder in the Linux file system using ``mkdir -p ~/tmp/sage`` and use it for building by ``eval SAGE_BUILD_DIR="~/tmp/sage" make``.
Also see the `related Github issue <https://github.com/pypa/packaging-problems/issues/258>`_ for other workarounds.

WSL post-installation notes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the installation is complete, you may be interested in :ref:`sec-launching-wsl-post-installation`.


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


Tcl/Tk (and system's Python)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to use `Tcl/Tk <https://www.tcl.tk/>`_ libraries in Sage, and you
are going to use your OS's Python3 as Sage's Python, you merely need to install
its **Tkinter** module.  On Linux systems, it is usually provided by the
**python3-tk** or a similarly named (e.g. **python3-tkinter**) package,
which can be installed using::

    $ sudo apt-get install python3-tk

or similar commands.

Tcl/Tk (and Sage's own Python)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to use `Tcl/Tk <https://www.tcl.tk/>`_ libraries in Sage,
and you are going to build Sage's Python from source, you need to install
these, and the corresponding headers.
On Linux systems, these are usually provided by the **tk** and **tk-dev**
(or **tk-devel**) packages which can be installed using::

    $ sudo apt-get install tk tk-dev

or similar commands.


Sage's Python will then automatically recognize your system's install of Tcl/Tk.
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

does not raise an :class:`ImportError`, then it worked.


.. _build-from-source-step-by-step:

Installation steps
------------------

.. hint:: 

  The following steps use the classical ``./configure && make`` build
  process. The modern Meson build system is also supported, see
  :ref:`build-source-meson`.

#. Follow the procedure in the file `README.md <https://github.com/sagemath/sage/#readme>`_
   in ``SAGE_ROOT``.

#. If you wish to prepare for having to build Sage in an environment
   without sufficient Internet connectivity:

   - After running ``configure``, you can use ``make download`` to force
     downloading packages before building. After this, the packages
     are in the subdirectory ``upstream``.

   - Alternatively, instead of cloning the git repository, you
     can download a self-contained release tarball for any
     stable release from the Sage project's
     `GitHub Releases <https://github.com/sagemath/sage/releases>`_.
     Use the file named ``sage-x.y.tar.gz`` (1.25 GB as of Sage 10.2)
     in the Release Assets, which contains a prepopulated subdirectory
     ``upstream``.

     After downloading the source tarball ``sage-x.y.tar.gz`` into
     a directory ``~/sage/``::

       $ cd ~/sage/
       $ tar xf sage-x.y.tar.gz  # adapt x.y; takes a while

     This creates the subdirectory ``sage-x.y``. Now change into it::

       $ cd sage-x.y/  # adapt x.y

     .. note::

        On Windows, it is crucial that you unpack the source tree from the
        WSL `bash` using the WSL `tar` utility and not using other
        Windows tools (including mingw).

        This is because the Sage source tree contains symbolic links, and the
        build will not work if Windows line endings rather than UNIX
        line endings are used.

   - The Sage mirrors also provide such self-contained tarballs
     for all `stable releases <https://www.sagemath.org/download-source.html>`_
     and additionally for all `development releases
     <https://www.sagemath.org/download-latest.html>`_.

#. Additional remarks:
   You do not need to be logged in as root, since no files are
   changed outside of the :file:`SAGE_ROOT` directory.
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
   (Linux, macOS, Solaris, OpenSolaris, or any other system),
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
     directory in your :envvar:`PATH`) to :sage_root:`sage`::

         $ ln -s /path/to/sage_root/sage /usr/local/bin/sage

     Now simply typing ``sage`` from any directory should be sufficient to run
     Sage.

   - Copy :sage_root:`sage` to a location in your :envvar:`PATH`.
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
     :sage_root:`sage`.
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
directory to build Sage and its documentation in HTML format, suitable for
viewing in a web browser.

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
(Note though that setting :envvar:`MAKEFLAGS`, as described below, can
significantly speed up the process.)
Building Sage involves building many packages, each of which has its own
compilation instructions.


Standard environment controlling the build process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here are some of the more commonly used variables affecting the build process:

.. envvar:: MAKEFLAGS

  This variable can be set to tell the ``make`` program to build things in
  parallel. Set it to ``-jNUM`` to run ``NUM`` jobs in parallel when building.
  Add ``-lNUM`` to tell make not to spawn more processes when the load exceeds
  ``NUM``.

  A good value for this variable is ``MAKEFLAGS="-j$(nproc) -l$(nproc).5"`` on
  Linux and ``MAKEFLAGS="-j$(sysctl -n hw.ncpu) -l$(sysctl -n hw.ncpu).5"`` on
  macOS. This instructs make to use all the execution threads of your CPU while
  bounding the load if there are other processes generating load. If your
  system does not have a lot of RAM, you might want to choose lower limits, if
  you have lots of RAM, it can sometimes be beneficial to set these limits
  slightly higher.

  Note that some parts of the SageMath build system do not respect this
  variable, e.g., when ninja gets invoked, it figures out the number of
  processes to use on its own so the number of processes and the system load
  you see might exceed the number configured here.

  See the manual page for GNU ``make``: `Command-line options
  <https://www.gnu.org/software/make/manual/make.html#Options-Summary>`_
  and `Parallel building
  <https://www.gnu.org/software/make/manual/make.html#Parallel>`_.

.. envvar:: V

  If set to ``0``, silence the build.  Instead of showing a detailed
  compilation log, only one line of output is shown at the beginning
  and at the end of the installation of each Sage package.  To see
  even less output, use::

    $ make -s V=0

  (Note that the above uses the syntax of setting a Makefile variable.)

.. envvar:: CC

  While some programs allow you to use this to specify your C
  compiler, **not every Sage package recognizes this**.
  If GCC is installed within Sage, :envvar:`CC` is ignored and Sage's ``gcc``
  is used instead.

.. envvar:: CPP

  Similarly, this will set the C preprocessor for some Sage
  packages, and similarly, using it is likely quite risky.
  If GCC is installed within Sage, :envvar:`CPP` is ignored and Sage's ``cpp``
  is used instead.

.. envvar:: CXX

  Similarly, this will set the C++ compiler for some Sage
  packages, and similarly, using it is likely quite risky.
  If GCC is installed within Sage, :envvar:`CXX` is ignored and Sage's ``g++``
  is used instead.

.. envvar:: FC

  Similarly, this will set the Fortran compiler.
  This is supported by all Sage packages which have Fortran code.
  However, for historical reasons, the value is hardcoded during the initial
  ``make`` and subsequent changes to ``$FC`` might be ignored (in which case,
  the original value will be used instead).
  If GCC is installed within Sage, :envvar:`FC` is ignored and Sage's
  ``gfortran`` is used instead.

.. envvar:: CFLAGS
.. envvar:: CXXFLAGS
.. envvar:: FCFLAGS

  The flags for
  the C compiler, the C++ compiler and the Fortran compiler, respectively.
  The same comments apply to these: setting them may cause problems, because
  they are not universally respected among the Sage packages. Note
  also that ``export CFLAGS=""`` does not have the same effect as
  ``unset CFLAGS``. The latter is preferable.

.. envvar:: CPPFLAGS
.. envvar:: LDFLAGS
.. envvar:: CXXFLAG64
.. envvar:: LDFLAG64
.. envvar:: LD

  Similar comments apply to these compiler and linker flags.


Sage-specific environment variables controlling the build process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. envvar:: SAGE_SERVER

  The Sage source tarball already includes the sources for all standard
  packages, that is, it allows you to build Sage without internet
  connection. The git repository, however, does not contain the source
  code for third-party packages. Instead, it will be downloaded as
  needed (note: you can run ``make download`` to force downloading
  packages before building).

  If :envvar:`SAGE_SERVER` is set, the specified Sage mirror is contacted
  first. Note that Sage will search the directory
  ``SAGE_SERVER/spkg/upstream`` for upstream tarballs.

  If downloading a file from there fails or :envvar:`SAGE_SERVER` is not set,
  files will be attempted to download from release assets of the
  Sage GitHub repository.

  If that fails too, the Sage mirror network is contacted to determine
  the nearest mirrors.

  This sequence of operations is defined by the files in the directory
  :sage_root:`.upstream.d`.

.. envvar:: SAGE_NUM_THREADS

  If set to a number, then when rebuilding with ``sage -b`` or
  parallel doctesting with ``sage -t -p 0``, use at most this many
  threads.

  If this is not set, then determine the number of threads using the value of
  the :envvar:`MAKE` (see above) or :envvar:`MAKEFLAGS` environment variables.
  If none of these specifies a number of jobs,

  - ``sage -b`` only uses one thread

  - ``sage -t -p 0`` uses a default of the number of CPU cores, with a
    maximum of 8 and a minimum of 2.

  When ``sage -t -p`` runs under the control of the GNU ``make``
  jobserver, then Sage will request as most this number of job slots.

.. envvar:: SAGE_CHECK

  If set to ``yes``, then during the build process,
  or when installing packages manually,
  run the test suite for each package which has one, and stop with an error
  if tests are failing.  If set to ``warn``, then only a warning is printed
  in this case.
  See also :envvar:`SAGE_CHECK_PACKAGES`.

.. envvar:: SAGE_CHECK_PACKAGES

  If :envvar:`SAGE_CHECK` is set to ``yes``,
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

.. envvar:: SAGE_INSTALL_GCC

  **Obsolete, do not use, to be removed**

.. envvar:: SAGE_INSTALL_CCACHE

  By default Sage doesn't install :ref:`ccache <spkg_ccache>`,
  however by setting ``SAGE_INSTALL_CCACHE=yes`` Sage will install ccache.
  Because the Sage distribution is quite large, the maximum cache is set to 4G.
  This can be changed by running ``sage -sh -c "ccache --max-size=SIZE"``,
  where ``SIZE`` is specified in gigabytes, megabytes, or kilobytes by
  appending a "G", "M", or "K".

  Sage does not include the sources for ccache since it is an optional package.
  Because of this, it is necessary to have an Internet connection while
  building ccache for Sage, so that Sage can pull down the necessary
  sources.

.. envvar:: SAGE_DEBUG

  Controls debugging support. There are three different possible values:

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

.. envvar:: SAGE_PROFILE

  Controls profiling support. If this is set
  to ``yes``, profiling support is enabled where possible. Note that
  Python-level profiling is always available; this option enables
  profiling in Cython modules.

.. envvar:: SAGE_BUILD_DIR

  The default behavior is to build each spkg in a
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

.. envvar:: SAGE_KEEP_BUILT_SPKGS

  The default behavior is to delete each
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

.. envvar:: SAGE_FAT_BINARY

  To build binaries that will run on the
  widest range of target CPUs set this variable to ``yes`` before
  building Sage or configure with ``--enable-fat-binary``.
  This does not make the binaries relocatable, it only
  avoids newer CPU instruction set extensions. For relocatable (=can
  be moved to a different directory) binaries, you must use
  https://github.com/sagemath/binary-pkg

.. envvar:: SAGE_SUDO

  Set this to ``sudo -E`` or to any other
  command prefix that is necessary to write into a installation
  hierarchy (:envvar:`SAGE_LOCAL`) owned by root or another user.
  Note that this command needs to preserve environment variable
  settings (plain ``sudo`` does not).

  Not all Sage packages currently support :envvar:`SAGE_SUDO`.

  Therefore this environment variable is most useful when a system
  administrator wishes to install an additional Sage package that
  supports :envvar:`SAGE_SUDO`, into a root-owned installation
  hierarchy (:envvar:`SAGE_LOCAL`).


Environment variables controlling the documentation build
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. envvar:: SAGE_DOCBUILD_OPTS

  The value of this variable is passed as an
  argument to ``sage --docbuild all html`` or ``sage --docbuild all pdf`` when
  you run ``make``, ``make doc``, or ``make doc-pdf``.  For example:

  - add ``--no-plot`` to this variable to avoid building the graphics coming from
    the ``.. PLOT`` directive within the documentation,

  - add ``--no-preparsed-examples`` to only show the original Sage code of
    "EXAMPLES" blocks, suppressing the tab with the preparsed, plain Python
    version, or

  - add ``--include-tests-blocks`` to include all "TESTS" blocks in the reference
    manual.

  Run ``sage --docbuild help`` to see the full list of options.

.. envvar:: SAGE_SPKG_INSTALL_DOCS

  If set to ``yes``, then install
  package-specific documentation to
  :file:`$SAGE_ROOT/local/share/doc/PACKAGE_NAME/` when an spkg is installed.
  This option may not be supported by all spkgs. Some spkgs might also assume
  that certain programs are available on the system (for example, ``latex`` or
  ``pdflatex``).

.. envvar:: SAGE_USE_CDNS

  If set to ``yes``, then build the documentation
  using CDNs (Content Distribution Networks) for scripts necessary for HTML
  documentation, such as `MathJax <https://www.mathjax.org/>`_.

.. envvar:: SAGE_LIVE_DOC

  If set to ``yes``, then build live Sage
  documentation. If the ``Make live`` button on any webpage of the live doc is
  clicked, every example code gets a `CodeMirror <https://codemirror.net>`_
  code cell runnable via `Thebe <https://thebe.readthedocs.io/en/stable/>`_.
  Thebe is responsible in sending the code to the Sage computing environment
  built by `Binder <https://mybinder.org/>`_ and showing the output result.
  The Sage computing environment can be specified to either a Binder repo or a
  local Jupyter server. The environment variable :envvar:`SAGE_JUPYTER_SERVER`
  is used for this purpose.

.. envvar:: SAGE_JUPYTER_SERVER

  Set this to either ``binder``, ``binder:repo`` with ``repo``
  specifying a Binder repo or the URL to a local Jupyter server.

  - ``binder`` refers to `Sage's official Binder repo
    <https://github.com/sagemath/sage-binder-env>`_. This is assumed if the
    environment variable :envvar:`SAGE_JUPYTER_SERVER` is not set.

  - ``binder:repo`` specifies a Binder repo with ``repo``, which is a GitHub
    repository name, optionally added with a branch name with ``/`` separator.

  - To use a local Jupyter server instead of Binder, then set the URL to
    :envvar:`SAGE_JUPYTER_SERVER` and the secret token to environment variable
    :envvar:`SAGE_JUPYTER_SERVER_TOKEN`, which can be left unset if the default
    token  ``secret`` is used. If the live doc was built with
    ``SAGE_JUPYTER_SERVER=http://localhost:8889``, run a local Jupyter server
    by

    .. CODE-BLOCK:: bash

        ./sage --notebook=jupyterlab \
               --ServerApp.token='secret' \
               --ServerApp.allow_origin='null' \
               --ServerApp.disable_check_xsrf=true \
               --ServerApp.port=8889 \
               --ServerApp.open_browser=false

    before opening the Sage documentation webpage.


Environment variables dealing with specific Sage packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. envvar:: SAGE_MATPLOTLIB_GUI

  If set to anything non-empty except ``no``,
  then Sage will attempt to build the graphical backend when it builds the
  matplotlib package.

.. envvar:: OPENBLAS_CONFIGURE

  Adds additional configuration flags for
  the OpenBLAS package that gets added to the ``make`` command. (see :issue:`23272`)

.. envvar:: PARI_CONFIGURE

  Use this to pass extra parameters to
  PARI's ``Configure`` script, for example to specify graphics
  support (which is disabled by default). See the file
  :file:`build/pkgs/pari/spkg-install.in` for more information.

.. envvar:: SAGE_TUNE_PARI

  If yes, enable PARI self-tuning. Note that
  this can be time-consuming. If you set this variable to "yes", you
  will also see this: ``WARNING: Tuning PARI/GP is unreliable. You may
  find your build of PARI fails, or PARI/GP does not work properly
  once built. We recommend to build this package with
  SAGE_CHECK="yes".``

.. envvar:: PARI_MAKEFLAGS

  The value of this variable is passed as an
  argument to the ``$MAKE`` command when compiling PARI.


Environment variables dealing with doctesting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. envvar:: SAGE_TIMEOUT

  Used for Sage's doctesting: the number of seconds
  to allow a doctest before timing it out.
  If this isn't set, the default is 300 seconds (5 minutes).

.. envvar:: SAGE_TIMEOUT_LONG

  Used for Sage's doctesting: the number of
  seconds to allow a doctest before timing it out, if tests are run using
  ``sage -t --long``.
  If this isn't set, the default is 1800 seconds (30 minutes).

.. envvar:: SAGE_TEST_GLOBAL_ITER
.. envvar:: SAGE_TEST_ITER

  These can
  be used instead of passing the flags ``--global-iterations`` and
  ``--file-iterations``, respectively, to ``sage -t``. Indeed, these
  variables are only used if the flags are unset. Run ``sage -t -h``
  for more information on the effects of these flags (and therefore
  these variables).


Environment variables set within Sage environments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

    - :envvar:`SAGE_TIMEOUT_VALGRIND` -- used for Sage's doctesting: the
      number of seconds to allow a doctest before timing it out, if tests
      are run using ``??``.  If this isn't set, the default is 1024*1024
      seconds.

    - :envvar:`SAGE_VALGRIND` -- trigger black magic in Python.

    - :envvar:`SAGE_MEMCHECK_FLAGS`, :envvar:`SAGE_MASSIF_FLAGS`,
      :envvar:`SAGE_CACHEGRIND_FLAGS`, :envvar:`SAGE_OMEGA_FLAGS` - flags
      used when using valgrind and one of the tools "memcheck", "massif",
      "cachegrind", or "omega"
    ***************************************************************************


Installation in a multiuser environment
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


Upgrading the system and upgrading Sage
---------------------------------------

Caveats when upgrading system packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When Sage has been installed from source, it will make use of various system
packages; in particular, it will link to shared libraries provided by
the system.

The system's package manager does not keep track of the applications that
make use of the shared libraries.  Therefore indiscriminate upgrades of
system packages can break a Sage installation.

This can always be fixed by a full rebuild::

  $ make distclean && make build

But this time-consuming step can often be avoided by just reinstalling a
few packages. The command ``make -j list-broken-packages`` assists with
this::

  $ make -j list-broken-packages
  make --no-print-directory auditwheel_or_delocate-no-deps
  ...
  # Checking .../local/var/lib/sage/installed/bliss-0.73+debian-1+sage-2016-08-02.p0
  ...
  Checking shared library file '.../local/lib/libumfpack.dylib'
  Checking shared library file '.../local/var/tmp/sage/build/suitesparse-5.10.1/src/lib/libsliplu.1.0.2.dylib'
  Error during installcheck of 'suitesparse': .../local/var/tmp/sage/build/suitesparse-5.10.1/src/lib/libsliplu.1.0.2.dylib
  ...
  Uninstall broken packages by typing:

      make lcalc-SAGE_LOCAL-uninstall;
      make ratpoints-SAGE_LOCAL-uninstall;
      make r-SAGE_LOCAL-uninstall;
      make suitesparse-SAGE_LOCAL-uninstall;

After running the suggested commands, run::

  $ make build


Upgrading Sage using a separate git worktree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When you have a working installation of Sage built from source and wish to
try out a new version, we strongly recommend to use a separate
`git worktree <https://git-scm.com/docs/git-worktree>`_, so that you
can keep using your existing installation when something goes wrong.

Start from the directory created when you used ``git clone``, perhaps
``~/sage/sage/``. Let's verify that this is indeed a git repository by
looking at the hidden ``.git`` subdirectory. It will looks like this,
but the exact contents can vary::

  [alice@localhost sage]$ ls .git
  COMMIT_EDITMSG HEAD           branches       description    gitk.cache
  index          logs           packed-refs    FETCH_HEAD     ORIG_HEAD
  config         hooks          info           objects        refs

Good. Now let's see what worktrees already exist::

  [alice@localhost sage]$ git worktree list
  /home/alice/sage/sage                     c0ffeefe10 [master]

We see just one line, the directory created when you used ``git clone``.
We will call this the "main worktree" from now on. Next to the directory,
you can see the abbreviated commit sha and the name of the branch that
we're on (``master``).

To try out a new version of Sage, let's fetch it first from the main
repository::

  [alice@localhost sage]$ git fetch upstream 10.3.beta8
  From https://github.com/sagemath/sage
   * tag                     10.3.beta8 -> FETCH_HEAD

Now let's create a new worktree. We need a name for it; it should
start with ``worktree-`` but can be anything after that. Experience
shows that worktrees are often repurposed later, and because a
directory containing a Sage installation cannot be moved without
breaking the installation in it, it may be a good idea to choose
a memorable name without much meaning::

  [alice@localhost sage]$ git worktree add worktree-purple FETCH_HEAD
  Preparing worktree (detached HEAD 30b3d78fac)
  Updating files: 100% (11191/11191), done.
  HEAD is now at 30b3d78fac Updated SageMath version to 10.3.beta8

We now have a subdirectory ``worktree-purple``. This is a
"linked worktree"::

  [alice@localhost sage]$ git worktree list
  /home/alice/sage/sage                     c0ffeefe10 [master]
  /home/alice/sage/sage/worktree-purple     30b3d78fac (detached HEAD)
  [alice@localhost sage]$ cd worktree-purple
  [alice@localhost worktree-purple]$ cat VERSION.txt
  SageMath version 10.3.beta8, Release Date: 2024-02-13

All worktrees created in this way share the same repository,
so they have access to all branches::

  [alice@localhost worktree-purple]$ git --no-pager branch -v
  * (no branch) 30b3d78fac Updated SageMath version to 10.3.beta8
  + master      2a9a4267f9 Updated SageMath version to 10.2

In fact, ``.git`` here is not a directory, just a hidden
file::

  [alice@localhost worktree-purple]$ ls -l .git
  -rw-r--r--  1 alice  staff  59 Feb 20 18:16 .git

In the new worktree, we now build Sage from scratch. This
is completely independent of and will not disrupt your
existing working installation in the main worktree.

We will refer again to the step-by-step instructions
from the file
`README.md <https://github.com/sagemath/sage/#readme>`_.
Our worktree ``worktree-purple`` is the ``SAGE_ROOT``
for this purpose.

One thing that we can share between worktrees without
worry is the directory ``upstream``, where Sage caches
downloaded archives of packages. To have the new worktree
share it with the main worktree, let's create a symbolic
link. This is an optional step that will avoid
re-downloading files that you already have::

  [alice@localhost worktree-purple]$ ln -s ../upstream/ .

Now let's build Sage, starting with the step::

  [alice@localhost worktree-purple]$ make configure

Refer to the file `README.md <https://github.com/sagemath/sage/#readme>`_
for the following steps.
