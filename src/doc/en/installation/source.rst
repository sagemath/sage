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

General requirements
~~~~~~~~~~~~~~~~~~~~

This section details the technical prerequisites needed on all platforms. See
also the `System-specific requirements`_ below.

Disk space and memory
^^^^^^^^^^^^^^^^^^^^^

Your computer comes with at least 6 GB of free disk space.
It is recommended to have at least 2 GB of RAM, but you might get away
with less (be sure to have some swap space in this case).

Command-line tools
^^^^^^^^^^^^^^^^^^

In addition to standard :wikipedia:`POSIX <POSIX>` utilities
and the :wikipedia:`bash <Bash_(Unix_shell)>` shell,
the following standard command-line development tools must be installed on your
computer:

- A **C/C++ compiler**: GCC versions 6.3 to 12.x are supported.
  Clang (LLVM) is also supported.
  See also `Using alternative compilers`_.
- **make**: GNU make, version 3.80 or later. Version 3.82 or later is recommended.
- **m4**: GNU m4 1.4.2 or later (non-GNU or older versions might also work).
- **perl**: version 5.8.0 or later.
- **ar** and **ranlib**: can be obtained as part of GNU binutils.
- **tar**: GNU tar version 1.17 or later, or BSD tar.
- **python**: Python 3.4 or later, or Python 2.7.
  (This range of versions is a minimal requirement for internal purposes of the SageMath
  build system, which is referred to as ``sage-bootstrap-python``.)

Other versions of these may work, but they are untested.


Fortran and compiler suites
###########################

Sage installation also needs a Fortran compiler.  It is determined
automatically whether Sage's GCC package, or just its part containing
Fortran compiler ``gfortran`` needs to be installed. This can be
overwritten by running ``./configure`` with option
``--without-system-gcc``.

Officially we support
gfortran from `GNU Compiler Collection (GCC) <https://gcc.gnu.org/>`_.
If C and C++ compilers also come from there (i.e., gcc and g++), their versions
should match.
Alternatively, one may use C and C++ compilers from
`Clang: a C language family frontend for LLVM <https://clang.llvm.org/>`_,
and thus  matching versions of
clang, clang++ , along with a recent gfortran. (Flang (or other LLVM-based
Fortran compilers) are not officially supported, however it is possible to
to build Sage using flang, with some extra efforts needed to set various flags;
this is work in progress at the moment (May 2019)).

Therefore, if you plan on using your own GCC compilers, then make sure that
their versions match.

To force using specific compilers, set environment variables ``CC``,
``CXX``, and ``FC`` (for C, C++, and Fortran compilers, respectively)
to the desired values, and run ``./configure``. For example,
``./configure CC=clang CXX=clang++ FC=gfortran`` will configure Sage
to be built with Clang C/C++ compilers and Fortran compiler
``gfortran``.

Alternatively, Sage includes a GCC package, so that C, C++ and Fortran
compilers will be built when the build system detects that it is needed,
e.g., non-GCC compilers, or
versions of the GCC compilers known to miscompile some components of Sage,
or simply a missing Fortran compiler.
In any case, you always need at least a C/C++ compiler to build the GCC
package and its prerequisites before the compilers it provides can be used.

Note that you can always override this behavior through the configure
options ``--without-system-gcc`` and ``--with-system-gcc``, see
:ref:`section_compilers`.

There are some known problems with old assemblers, in particular when
building the ``ecm`` and ``fflas_ffpack`` packages. You should ensure
that your assembler understands all instructions for your
processor. On Linux, this means you need a recent version of
``binutils``; on macOS you need a recent version of Xcode.

Python for venv
^^^^^^^^^^^^^^^

By default, Sage will try to use system's `python3` to set up a virtual
environment, a.k.a. `venv <https://docs.python.org/3.10/library/venv.html>`_
rather than building a Python 3 installation from scratch.
Use the ``configure`` option ``--without-system-python3`` in case you want Python 3
built from scratch.

Sage will accept versions 3.8.x to 3.10.x.

You can also use ``--with-python=/path/to/python3_binary`` to tell Sage to use
``/path/to/python3_binary`` to set up the venv. Note that setting up venv requires
a number of Python modules to be available within the Python in question. Currently,
for Sage 9.6, these modules are as follows: ``sqlite3``, ``ctypes``, ``math``,
``hashlib``, ``crypt``, ``socket``, ``zlib``, ``distutils.core``, ``ssl`` -
they will be checked for by the ``configure`` script.

Other notes
^^^^^^^^^^^

After extracting the Sage source tarball, the subdirectory :file:`upstream`
contains the source distributions for everything on which Sage depends.

If cloned from a git repository, the upstream tarballs will be downloaded,
verified, and cached as part of the Sage installation process.
We emphasize that all of this software is included with Sage, so you do not
have to worry about trying to download and install any one of these packages
(such as Python, for example) yourself.

When the Sage installation program is run,
it will check that you have each of the above-listed prerequisites,
and inform you of any that are missing, or have unsuitable versions.

System-specific requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On macOS, there are various developer tools needed which may require
some registration on Apple's developer site; see
:ref:`section_macprereqs`.

On Redhat-derived systems not all perl components are installed by
default and you might have to install the ``perl-ExtUtils-MakeMaker``
package.

On Cygwin, the ``lapack`` and ``liblapack-devel`` packages are required.

Installing prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~

To check if you have the above prerequisites installed, for example ``perl``,
type::

    $ command -v perl

or::

    $ which perl

on the command line. If it gives an error (or returns nothing), then
either ``perl`` is not installed, or it is installed but not in your
:wikipedia:`PATH <PATH_%28variable%29>`.

.. _sec-installation-from-sources-linux-recommended-installation:

Linux recommended installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On Linux systems (e.g., Ubuntu, Redhat, etc), ``ar`` and ``ranlib`` are in the
`binutils <https://www.gnu.org/software/binutils/>`_ package.
The other programs are usually located in packages with their respective names.
Assuming you have sufficient privileges, you can install the ``binutils`` and
other necessary/standard components. The lists provided below are longer than
the minimal prerequisites, which are basically ``binutils``, ``gcc``/``clang``, ``make``,
``tar``, but there is no real need to build compilers and other standard tools
and libraries on a modern Linux system, in order to be able to build Sage.
If you do not have the privileges to do this, ask your system administrator to
do this, or build the components from source code.
The method of installing additional software varies from distribution to
distribution, but on a `Debian <https://www.debian.org/>`_ based system (e.g.
`Ubuntu <https://www.ubuntu.com/>`_ or `Mint <https://www.linuxmint.com/>`_),
you would use
:wikipedia:`apt-get <Advanced_Packaging_Tool>`.

On Debian ("buster" or newer) or Ubuntu ("bionic" or newer):

.. literalinclude:: debian.txt

On Fedora / Redhat / CentOS:

.. literalinclude:: fedora.txt

On Arch Linux:

.. literalinclude:: arch.txt

In addition to these, if you don't want Sage to build optional packages that might
be available from your OS, cf. the growing list of such packages on :trac:`27330`,
install on Debian ("buster" or newer) or Ubuntu ("bionic" or newer):

.. literalinclude:: debian-optional.txt

On Fedora / Redhat / CentOS:

.. literalinclude:: fedora-optional.txt

On Arch Linux:

.. literalinclude:: arch-optional.txt

On other Linux systems, you might use
:wikipedia:`rpm <RPM_Package_Manager>`,
:wikipedia:`yum <Yellowdog_Updater,_Modified>`,
or other package managers.

.. _section_macprereqs:

macOS prerequisite installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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



macOS recommended installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although Sage can in theory build its own version of gfortran, this
can take a while, and the process fails on some recent versions of
OS X. So instead you can install your own copy. One advantage of this
is that you can install it once, and it will get used every time you
build Sage, rather than building gfortran every time.

One way to do that is with the `Homebrew package manager
<https://brew.sh>`_. Install Homebrew as their web page describes, and
then the command ::

    $ brew install gcc

will install Homebrew's gcc package, which includes gfortran. Sage
will also use other Homebrew packages, if they are present. You can
install the following:

.. literalinclude:: homebrew.txt

Some Homebrew packages are installed "keg-only," meaning that they are
not available in standard paths. To make them accessible when building
Sage, run ::

    $ source SAGE_ROOT/.homebrew-build-env

(replacing ``SAGE_ROOT`` by Sage's home directory). You can add a
command like this to your shell profile if you want the settings to
persist between shell sessions.

Some additional optional packages are taken care of by:

.. literalinclude:: homebrew-optional.txt


.. _section_cygwinprereqs:

Cygwin prerequisite installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sage can be built only on the 64-bit version of Cygwin.  See
the file `README.md <https://github.com/sagemath/sage/#readme>`_
in ``SAGE_ROOT`` for the most up-to-date instructions for building Sage
on Cygwin.

Although it is possible to install Sage's dependencies using the Cygwin
graphical installer, it is recommended to install the `apt-cyg
<https://github.com/transcode-open/apt-cyg>`_ command-line package
installer, which is used for the remainder of these instructions.  To
run ``apt-cyg``, you must have already installed (using the graphical
installer) the following packages at a minimum::

    bzip2 coreutils gawk gzip tar wget

With the exception of ``wget`` most of these are included in the default
package selection when you install Cygwin.  Then, to install ``apt-cyg``
run::

    $ curl -OL https://rawgit.com/transcode-open/apt-cyg/master/apt-cyg
    $ install apt-cyg /usr/local/bin
    $ rm -f apt-cyg

To install the current set of system packages known to work for building
Sage, run:

.. literalinclude:: cygwin.txt

Optional packages that are also known to be installable via system packages
include:

.. literalinclude:: cygwin-optional.txt

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

Specific notes for ``make`` and ``tar``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On macOS, the system-wide BSD ``tar`` supplied will build Sage, so there is no
need to install the GNU ``tar``.

.. _section_compilers:

Using alternative compilers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sage developers tend to use fairly recent versions of GCC.
Nonetheless, the Sage build process on Linux
should succeed with any reasonable C/C++ compiler;
(we do not recommend GCC older than version 5.1).
This is because Sage will build GCC first (if needed) and then use that newly
built GCC to compile Sage.

If you don't want this and want to try building Sage with a different set of
compilers,
you need to pass Sage's ``./configure`` compiler names, via environment
variables ``CC``, ``CXX``, and ``FC``, for C, C++, and Fortran compilers,
respectively, e.g. if you C compiler is ``clang``, your C++ compiler is ``clang++``,
and your Fortran compiler is ``flang`` then you would need to run::

    $ CC=clang CXX=clang++ FC=flang ./configure

before running ``make``. It is recommended that you inspect the output of ``./configure``
in order to check that Sage will not try to build GCC. Namely, there should be lines like::

       gcc-7.2.0 will not be installed (configure check)
       ...
       gfortran-7.2.0 will not be installed (configure check)

indicating that Sage will not attempt to build ``gcc/g++/gfortran``.

If you are interested in working on support for commercial compilers from
`HP <http://docs.hp.com/en/5966-9844/ch01s03.html>`_,
`IBM <http://www-01.ibm.com/software/awdtools/xlcpp/>`_,
`Intel <http://software.intel.com/en-us/articles/intel-compilers/>`_,
`Sun/Oracle <http://www.oracle.com/technetwork/server-storage/solarisstudio/overview/index.html>`_,
etc,
please email the sage-devel mailing list at https://groups.google.com/group/sage-devel.


Additional software
-------------------

Recommended programs
~~~~~~~~~~~~~~~~~~~~

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
~~~~~~

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
~~~~~~~~~~~~~~~~~

#. Follow the procedure in the file `README.md <https://github.com/sagemath/sage/#readme>`_
   in ``SAGE_ROOT``.

#. Additional remarks:
   You do not need to be logged in as root, since no files are
   changed outside of the :file:`sage-x.y` directory.
   In fact, **it is inadvisable to build Sage as root**, as the root account
   should only be used when absolutely necessary and mistyped commands can have
   serious consequences if you are logged in as root.

   Typing ``make`` performs the usual steps for each Sage's dependency,
   but installs all the resulting files into the local build tree.
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
   Install optional Sage packages and databases.
   Type ``sage --optional`` to see a list of them (this requires an Internet
   connection), or visit https://www.sagemath.org/packages/optional/.
   Then type ``sage -i <package-name>`` to automatically download and install
   a given package.

#. Have fun! Discover some amazing conjectures!

Rebasing issues on Cygwin
~~~~~~~~~~~~~~~~~~~~~~~~~

Building on Cygwin will occasionally require "rebasing" ``dll`` files.
Sage provides some scripts, located in :file:`$SAGE_LOCAL/bin`, to do so:

- ``sage-rebaseall.sh``, a shell script which calls Cygwin's ``rebaseall``
  program.
  It must be run within a ``dash`` shell from the :envvar:`SAGE_ROOT` directory
  after all other Cygwin processes have been shut down and needs write-access
  to the system-wide rebase database located at :file:`/etc/rebase.db.i386`,
  which usually means administrator privileges.
  It updates the system-wide database and adds Sage dlls to it, so that
  subsequent calls to ``rebaseall`` will take them into account.
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

Finally, on Cygwin, one should also avoid the following:

- building in home directories of Windows domain users;
- building in paths with capital letters
  (see :trac:`13343`, although there has been some success doing so).


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

Sage uses the following environment variables when it runs:

- :envvar:`DOT_SAGE` - this is the directory, to which the user has read and
  write access, where Sage stores a number of files.
  The default location is :file:`$HOME/.sage/`.

- :envvar:`SAGE_STARTUP_FILE` - a file including commands to be executed every
  time Sage starts.
  The default value is :file:`$DOT_SAGE/init.sage`.

- :envvar:`BROWSER` - on most platforms, Sage will detect the command to
  run a web browser, but if this doesn't seem to work on your machine, set this
  variable to the appropriate command.

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


**This page was last updated in May 2022 (Sage 9.7).**
