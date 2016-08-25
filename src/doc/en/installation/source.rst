.. comment:
    ***************************************************************************
    If you alter this document, please change the last line:
    **This page was last updated in MONTH YEAR (Sage X.Y).**
    ***************************************************************************

.. _sec-installation-from-sources:

Install from Source Code
========================

More familiarity with computers may be required to build Sage from
the `source code <http://en.wikipedia.org/wiki/Source_code>`_.
If you do have all the :ref:`pre-requisite tools <section-prereqs>`,
the process should be completely
painless, basically consisting in extracting the source tarball and typing
``make``.  It can take your computer a while to build Sage from the source code,
although the procedure is fully automated and should need no human
intervention.

Building Sage from the source code has the major advantage that your install
will be optimised for your particular computer and should therefore offer
better performance and compatibility than a binary install.
Moreover, it offers you full development capabilities:
you can change absolutely any part of Sage or the programs on which it depends,
and recompile the modified parts.

`Download the Sage source code <http://www.sagemath.org/download-source.html>`_
or `check it out with git <https://github.com/sagemath/sage>`_ (see also.
`the developers guide <http://doc.sagemath.org/html/en/developer/manual_git.html#section-git-checkout>`_).
If you changed your mind, you can also download a
`binary distribution <http://www.sagemath.org/download.html>`_
for some operating systems.

.. contents:: Table of contents
   :depth: 2

Supported platforms
-------------------

See http://wiki.sagemath.org/SupportedPlatforms for the full list of platforms
on which Sage is supported and the level of support for these systems.

Sage is supported on a number of `Linux <http://en.wikipedia.org/wiki/Linux>`_,
Mac `OS X <http://www.apple.com/macosx/>`_ ,
Sun/Oracle `Solaris <http://www.oracle.com/solaris>`_ releases,
but not necessarily all versions of these operating systems.
There is no native version of Sage which installs on
`Microsoft Windows <http://en.wikipedia.org/wiki/Microsoft_Windows>`_,
although Sage can be used on Windows with the aid of a
`virtual machine <http://en.wikipedia.org/wiki/Virtual_machine>`_
or the `Cygwin <http://cygwin.com/>`_ Linux API layer.

On the `list of supported platforms <http://wiki.sagemath.org/SupportedPlatforms>`_,
you can find details about
`ports <http://en.wikipedia.org/wiki/Computer_port_%28software%29>`_
to other operating systems or processors which may be taking place.

.. _section-prereqs:

Prerequisites
-------------

General requirements
~~~~~~~~~~~~~~~~~~~~

This section details the technical prerequisites needed on all platforms. See
also the `System-specific requirements`_ below.

Disk space and memory
^^^^^^^^^^^^^^^^^^^^^

Your computer comes with at least 6 GB of free disk space running one of the
supported versions of an operating system listed at
http://wiki.sagemath.org/SupportedPlatforms.
It is recommended to have at least 2 GB of RAM, but you might get away
with less (be sure to have some swap space in this case).

Command-line tools
^^^^^^^^^^^^^^^^^^

In addition to standard `POSIX <http://en.wikipedia.org/wiki/POSIX>`_ utilities
and the `bash <http://en.wikipedia.org/wiki/Bash_(Unix_shell)>`_ shell,
the following standard command-line development tools must be installed on your
computer:

- A **C/C++ compiler**: Since Sage builds its own GCC if needed,
  a wide variety of C/C++ compilers is supported.
  Many GCC versions work,
  from as old as version 3.4.3 to the most recent release.
  Clang also works.
  On Solaris systems, the Sun compiler should also work.
  See also `Using alternative compilers`_.
- **make**: GNU make, version 3.80 or later. Version 3.82 or later is recommended.
- **m4**: GNU m4 1.4.2 or later (non-GNU or older versions might also work).
- **perl**: version 5.8.0 or later.
- **ar** and **ranlib**: can be obtained as part of GNU binutils.
- **tar**: GNU tar version 1.17 or later, or BSD tar.
- **python**: Python >= 2.6.

Fortran and compiler suites
###########################

Sage also needs a Fortran compiler.
The only configuration currently supported is matching versions of the
C, C++ and Fortran compilers from the
`GNU Compiler Collection (GCC) <http://gcc.gnu.org/>`_.
Therefore, if you plan on using your own GCC compilers, then make sure that
their versions match.

Alternatively, Sage includes a GCC package, so that C, C++ and Fortran
compilers will be built when the build system detects that it is needed,
e.g., non-GCC compilers, or
versions of the GCC compilers known to miscompile some components of Sage,
or simply a missing Fortran compiler.
In any case, you always need at least a C/C++ compiler to build the GCC
package and its prerequisites before the compilers it provides can be used.

Note that you can always override this behavior through the environment
variable :envvar:`SAGE_INSTALL_GCC`, see :ref:`section_compilers` and
:ref:`section_envvar`.

Other notes
^^^^^^^^^^^

Although some of Sage is written in `Python <http://www.python.org/>`_, you do
not need Python pre-installed on your computer, since the Sage installation
includes virtually everything you need.

After extracting the Sage tarball, the subdirectory :file:`upstream`
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

On Mac OS X, there are various developer tools needed which may require
some registration on Apple's developer site; see
:ref:`section_macprereqs`.

On Redhat-derived systems not all perl components are installed by
default and you might have to install the **perl-ExtUtils-MakeMaker**
package.

On Cygwin, the **lapack** and **liblapack-devel** packages are required to
provide ATLAS support as the Sage package for ATLAS is not built by default.

Installing prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~

To check if you have the above prerequisites installed, for example ``perl``,
type::

    command -v perl

or::

    which perl

on the command line. If it gives an error (or returns nothing), then
either ``perl`` is not installed, or it is installed but not in your
`PATH <http://en.wikipedia.org/wiki/PATH_%28variable%29>`_.

Linux prerequisite installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On Linux systems (e.g., Ubuntu, Redhat, etc), ``ar`` and ``ranlib`` are in the
`binutils <http://www.gnu.org/software/binutils/>`_ package.
The other programs are usually located in packages with their respective names.
Assuming you have sufficient privileges, you can install the ``binutils`` and
other necessary components.
If you do not have the privileges to do this, ask your system administrator to
do this, or build the components from source code.
The method of installing additional software varies from distribution to
distribution, but on a `Debian <http://www.debian.org/>`_ based system (e.g.
`Ubuntu <http://www.ubuntu.com/>`_ or `Mint <http://www.linuxmint.com/>`_),
you would use
`apt-get <http://en.wikipedia.org/wiki/Advanced_Packaging_Tool>`_::

     # debian
     sudo apt-get install binutils gcc make m4 perl tar git

     # redhat
     sudo yum install binutils gcc make m4 perl tar git perl-ExtUtils-MakeMaker
     
to install all general requirements, or, if you don't want Sage to build its
own GCC::

     # debian
     sudo apt-get install binutils gcc g++ gfortran make m4 perl tar git

     # redhat
     sudo yum install binutils gcc gcc-c++ gcc-gfortran make m4 perl tar git perl-ExtUtils-MakeMaker
     
(This was tested on Ubuntu 12.04.2.)
On other Linux systems, you might use
`rpm <http://en.wikipedia.org/wiki/RPM_Package_Manager>`_,
`yum <http://en.wikipedia.org/wiki/Yellowdog_Updater,_Modified>`_,
or other package managers.

.. _section_macprereqs:

Mac OS X prerequisite installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On OS X systems, you need a recent version of
`Command Line Tools <http://developer.apple.com/downloads/index.action?=command%20line%20tools>`_.
It provides all the above requirements.

If you have already installed `Xcode <http://developer.apple.com/xcode/>`_
(which at the time of writing is freely available in the Mac App Store,
or through http://developer.apple.com/downloads/ provided you registered for an
Apple Developer account), you can install the command line tools from
there as well.

- With OS X Mavericks or Yosemite, run the command
  ``xcode-select --install`` from a Terminal window and click "Install"
  in the pop-up dialog box.

- Using OS X Mountain Lion or earlier, run Xcode, open its "Downloads"
  preference pane and install the command line tools from there.

- On pre-Lion OS X systems, the command line tools are not available as a
  separate download and you have to install the full-blown Xcode supporting your
  system version.

If you have not installed `Xcode <http://developer.apple.com/xcode/>`_
you can get these tools as a relatively small download, but it does require
a registration.

- First, you will need to register as an Apple Developer at
  http://developer.apple.com/register/.

- Having done so, you should be able to download it for free at
  http://developer.apple.com/downloads/index.action?=command%20line%20tools

- Alternately, https://developer.apple.com/opensource/ should have a link
  to Command Line Tools.

Other platforms
^^^^^^^^^^^^^^^

On Solaris, you would use ``pkgadd`` and on OpenSolaris ``ipf`` to install
the necessary software.

On Cygwin, you would use the ``setup.exe`` program.
As on Linux systems, ``ar`` and ``ranlib`` are provided by the ``binutils`` package.
As far as compilers are concerned, you should either install matching versions
of the ``gcc4-core``, ``gcc4-g++``, and ``gcc4-gfortran`` packages, or
the ``gcc4-core`` package alone if you plan on using Sage's own GCC.

On other systems, check the documentation for your particular operating system.

Specific notes for ``make`` and ``tar``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On OS X, the system-wide BSD ``tar`` supplied will build Sage, so there is no
need to install the GNU ``tar``.

On Solaris or OpenSolaris, the Sun/Oracle versions of ``make`` and ``tar`` are
unsuitable for building Sage.
Therefore, you must have the GNU versions of ``make`` and ``tar`` installed and
they must be the first ``make`` and ``tar`` in your :envvar:`PATH`.

On Solaris 10, a version of GNU ``make`` may be found at
:file:`/usr/sfw/bin/gmake`,
but you will need to copy it somewhere else and rename it to ``make``.
The same is true for GNU ``tar``; a version of GNU ``tar`` may be found at
:file:`/usr/sfw/bin/gtar`,
but it will need to be copied somewhere else and renamed to ``tar``.
It is recommended to create a directory :file:`$HOME/bins-for-sage` and to put
the GNU versions of ``tar`` and ``make`` in that directory.
Then ensure that :file:`$HOME/bins-for-sage` is first in your :envvar:`PATH`.
That's because Sage also needs :file:`/usr/ccs/bin` in your :envvar:`PATH` to
execute programs like ``ar`` and ``ranlib``, but :file:`/usr/ccs/bin` has the
Sun/Oracle versions of ``make`` and ``tar`` in it.

If you attempt to build Sage on AIX or HP-UX, you will need to install both
GNU ``make`` and GNU ``tar``.

.. _section_compilers:

Using alternative compilers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sage developers tend to use fairly recent versions of GCC.
Nonetheless, the Sage build process should succeed with any reasonable C/C++ compiler.
This is because Sage will build GCC first (if needed) and then use that newly
built GCC to compile Sage.

If you don't want this and want to try building Sage with a different set of
compilers,
you need to set the environment variable :envvar:`SAGE_INSTALL_GCC` to ``no``.
Make sure you have C, C++, and Fortran compilers installed!

Building all of Sage with Clang is currently not supported, see :trac:`12426`.

If you are interested in working on support for commerical compilers from
`HP <http://docs.hp.com/en/5966-9844/ch01s03.html>`_,
`IBM <http://www-01.ibm.com/software/awdtools/xlcpp/>`_,
`Intel <http://software.intel.com/en-us/articles/intel-compilers/>`_,
`Sun/Oracle <http://www.oracle.com/technetwork/server-storage/solarisstudio/overview/index.html>`_,
etc,
please email the sage-devel mailing list at http://groups.google.com/group/sage-devel.


Additional software
-------------------

Recommended programs
~~~~~~~~~~~~~~~~~~~~

The following programs are recommended.
They are not strictly required at build time or at run time,
but provide additional capablities:

- **dvipng**.
- **ffmpeg**.
- **ImageMagick**.
- **latex**: highly recommended.

It is highly recommended that you have
`Latex <http://en.wikipedia.org/wiki/LaTeX>`_
installed, but it is not required.

The most popular packaging is `TeX Live <www.tug.org/texlive/>`_ ,
which you can install locally inside Sage with the commands::

    sage -sh -c '$SAGE_ROOT/src/ext/texlive/texlive-install'

On Linux systems you can alternatively install your distribution's
texlive packages::

    sudo apt-get install texlive       # debian
    sudo yum install texlive           # redhat

or similar commands. In addition to the base texlive install you will
probably need a number of optional texlive packages, for example
country-specific babel packages for the localized Sage
documentation. The required texlive packages are listed in
``SAGE_ROOT/src/ext/texlive/package-list.txt``.

If you don't have either ImageMagick or ffmpeg, you won't be able to
view animations.
ffmpeg can produce animations in more different formats than ImageMagick,
and seems to be faster than ImageMagick when creating animated GIFs.
Either ImageMagick or dvipng is used for displaying some LaTeX output in the
Sage notebook.

Notebook additional features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the Sage notebook uses the
`HTTP <http://en.wikipedia.org/wiki/HTTP>`_
protocol when you type the command ``notebook()``.
To run the notebook in secure mode by typing ``notebook(secure=True)`` which
uses the `HTTPS <http://en.wikipedia.org/wiki/HTTPS>`_ protocol,
or to use `OpenID <http://en.wikipedia.org/wiki/OpenID>`_ authentication,
you need to follow specific installation steps described in
:ref:`section_notebook_ssl`.

Although all necessary components are provided through Sage optional packages,
i.e. you can install a local version of `OpenSSL <http://www.openssl.org>`_
by using Sage's **openssl** package and running ``sage -i openssl`` as suggested
in :ref:`section_notebook_ssl` (this requires an Internet connection),
you might prefer to install OpenSSL and the OpenSSL development headers
globally on your system.

On Linux systems, those are usually provided by the **libssl** and
**libssl-dev** packages and can be installed using::

    sudo apt-get install libssl libssl-dev

or similar commands.

Finally, if you intend to distribute the notebook load onto several Sage
servers, you will surely want to setup an
`SSH <http://en.wikipedia.org/wiki/SSH>`_ server and generate SSH keys.
This can be achieved using `OpenSSH <http://www.openssh.org>`_.

On Linux systems, the OpenSSH server, client and utilities are usually provided
by the **openssh-server** and **openssh-client** packages and can be installed
using::

    sudo apt-get install openssh-server openssh-client

or similar commands.

Tcl/Tk
~~~~~~

If you want to use `Tcl/Tk <http://www.tcl.tk/>`_ libraries in Sage,
you need to install the Tcl/Tk and its development headers before building
Sage.
Sage's Python will then automatically recognize your system's install of
Tcl/Tk.

On Linux systems, these are usually provided by the **tk** and **tk-dev**
(or **tk-devel**) packages which can be installed using::

    sudo apt-get install tk tk-dev

or similar commands.

If you installed Sage first, all is not lost. You just need to rebuild
Sage's Python and any part of Sage relying on it::

    sage -f python2  # rebuild Python
    make             # rebuild components of Sage depending on Python

after installing the Tcl/Tk development libraries as above.

If

.. skip

::

   sage: import _tkinter
   sage: import Tkinter

does not raise an ``ImportError``, then it worked.

.. _build-from-source-step-by-step:

Step-by-step installation procedure
-----------------------------------

General procedure
~~~~~~~~~~~~~~~~~

Installation from source is (potentially) very easy, because the distribution
contains (essentially) everything on which Sage depends.

Make sure there are **no spaces** in the path name for the directory in which
you build:
several of Sage's components will not build if there are spaces in the path.
Running Sage from a directory with spaces in its name will also fail.

#. Go to http://www.sagemath.org/download-source.html, select a mirror,
   and download the file :file:`sage-x.y.tar`.

   This tarfile contains the source code for Sage and the source for all
   programs on which Sage depends.
   Note that this file is not compressed; it's just a plain tarball (which
   happens to be full of compressed files).

   Download it into any directory you have write access to, preferably on a
   fast filesystem, avoiding
   `NFS <http://en.wikipedia.org/wiki/Network_File_System>`_ and the like.
   On personal computers, any subdirectory of your :envvar:`HOME` directory
   should do.
   The directory where you built Sage is **NOT** hardcoded.
   You should be able to safely move or rename that directory.
   (It's a bug if this is not the case.)

#. Extract the tarfile::

       tar xvf sage-x.y.tar

   This creates a directory :file:`sage-x.y`.

#. Change into that directory::

       cd sage-x.y

   This is Sage's home directory.
   It is also referred to as :envvar:`SAGE_ROOT` or the top level Sage
   directory.

#. Optional, but highly recommended:
   Read the :file:`README.txt` file there.

#. On OSX 10.4, OS 10.5, Solaris 10 and OpenSolaris, if you wish to build a
   64-bit version of Sage, assuming your computer and operating system are
   64-bit, type::

       export SAGE64=yes

   It should be noted that as of April 2011, 64-bit builds of Sage on both
   Solaris 10 and OpenSolaris are not very stable, so you are advised not to
   set :envvar:`SAGE64` to ``yes``.
   This will then create stable 32-bit versions of Sage.
   See http://wiki.sagemath.org/solaris for the latest information.

#. Start the build process::

       make

   or if your system supports multiprocessing and you want to use several
   processes to build Sage::

       MAKE='make -jNUM' make

   to tell the ``make`` program to run ``NUM`` jobs in parallel when building
   Sage. This compiles Sage and all its dependencies.

   .. NOTE::

      Mac OS X allows changing directories without using exact capitalization.
      Beware of this convenience when compiling for OS X. Ignoring exact
      capitalization when changing into :envvar:`SAGE_ROOT` can lead to build
      errors for dependencies requiring exact capitalization in path names.

   Note that you do not need to be logged in as root, since no files are
   changed outside of the :file:`sage-x.y` directory.
   In fact, **it is inadvisable to build Sage as root**, as the root account
   should only be used when absolutely necessary and mistyped commands can have
   serious consequences if you are logged in as root.
   There has been a bug reported (see :trac:`9551`) in Sage which would have
   overwritten a system file had the user been logged in as root.

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
   newsgroup at http://groups.google.com/group/sage-support.
   If the log files are very large (and many are), then don't paste the whole
   file, but make sure to include any error messages.
   It would also be helpful to include the type of operating system
   (Linux, OS X, Solaris, OpenSolaris, Cygwin, or any other system),
   the version and release date of that operating system and the version of
   the copy of Sage you are using.
   (There are no formal requirements for bug reports -- just send them;
   we appreciate everything.)

   See :ref:`section_make` for some targets for the ``make`` command,
   :ref:`section_envvar` for additional informatio on useful environment
   variables used by Sage,
   and :ref:`section_notebook_ssl` for additional instruction on how to build
   the notebook with SSL support.

#. To start Sage, you can now simply type from Sage's home directory::

       ./sage

   You should see the Sage prompt, which will look something like this::

       $ sage
       ----------------------------------------------------------------------
       | Sage Version 5.8, Release Date: 2013-03-15                         |
       | Type "notebook()" for the browser-based notebook interface.        |
       | Type "help()" for help.                                            |
       ----------------------------------------------------------------------
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
   report the problem, e.g., at http://groups.google.com/group/sage-support.

   After Sage has started, try a simple command::

       sage: 2 + 2
       4

   Or something slightly more complicated::

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
   name: `Mathematica <http://www.wolfram.com/mathematica/>`_ is invoked by
   calling ``math``, `Maple <http://www.maplesoft.com/>`_ by calling ``maple``,
   etc.
   The easiest way to change this name or perform other customizations is
   to create a redirection script in :file:`$SAGE_ROOT/local/bin`.
   Sage inserts this directory at the front of your :envvar:`PATH`, so your
   script may need to use an absolute path to avoid calling itself; also, your
   script should pass along all of its arguments.
   For example, a ``maple`` script might look like::

       #!/bin/sh

       exec /etc/maple10.2/maple.tty "$@"

#. Optional:
   There are different possibilities to make using Sage a little easier:

   - Make a symbolic link from :file:`/usr/local/bin/sage` (or another
     directory in your :envvar:`PATH`) to :file:`$SAGE_ROOT/sage`::

         ln -s /path/to/sage-x.y/sage /usr/local/bin/sage

     Now simply typing ``sage`` from any directory should be sufficient to run
     Sage.

   - Copy :file:`$SAGE_ROOT/sage` to a location in your :envvar:`PATH`.
     If you do this, make sure you edit the line::

         #SAGE_ROOT=/path/to/sage-version

     at the beginning of the copied ``sage`` script according to the direction
     given there to something like::

         SAGE_ROOT=<SAGE_ROOT>

     (note that you have to change ``<SAGE_ROOT>`` above!).
     It is best to edit only the copy, not the original.

   - For `KDE <http://www.kde.org/>`_ users, create a bash script called
     ``sage`` containing the lines
     (note that you have to change ``<SAGE_ROOT>`` below!)::

         #!/bin/bash

         konsole -T "sage" -e <SAGE_ROOT>/sage

     make it executable::

         chmod a+x sage

     and put it somewhere in your :envvar:`PATH`.

     You can also make a KDE desktop icon with this line as the command
     (under the Application tab of the Properties of the icon, which you get my
     right clicking the mouse on the icon).

   - On Linux and OS X systems, you can make an alias to
     :file:`$SAGE_ROOT/sage`.
     For example, put something similar to the following line in your
     :file:`.bashrc` file::

         alias sage=<SAGE_ROOT>/sage

     (Note that you have to change ``<SAGE_ROOT>`` above!)
     Having done so, quit your terminal emulator and restart it.
     Now typing ``sage`` within your terminal emulator should start Sage.

#. Optional:
   Install optional Sage packages and databases.
   Type ``sage --optional`` to see a list of them (this requires an Internet
   connection), or visit http://www.sagemath.org/packages/optional/.
   Then type ``sage -i <package-name>`` to automatically download and install
   a given package.

#. Optional:
   Run the ``install_scripts`` command from within Sage to create GAP, GP,
   Maxima, Singular, etc., scripts in your :envvar:`PATH`.
   Type ``install_scripts?`` in Sage for details.

#. Have fun! Discover some amazing conjectures!

.. _section_notebook_ssl:

Building the notebook with SSL support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Read this section if you are intending to run a Sage notebook server for
multiple users.

For security, you may wish users to access the server using the HTTPS protocol
(i.e., to run ``notebook(secure=True)``).
You also may want to use OpenID for user authentication.
The first of these requires you to install
`pyOpenSSL <http://pyopenssl.sourceforge.net/>`_,
and they both require OpenSSL.

If you have OpenSSL and the OpenSSL development headers installed on your
system, you can install pyOpenSSL by building Sage and then typing::

    ./sage -i pyopenssl

Alternatively, ``make ssl`` builds Sage and installs pyOpenSSL at once.
Note that these commands require Internet access.

If you are missing either OpenSSL or OpenSSL's development headers,
you can install a local copy of both into your Sage installation first.
Ideally, this should be done before installing Sage; otherwise, you should at
least rebuild Sage's Python, and ideally any part of Sage relying on it.
The procedure is as follows (again, with a computer connected to the
Internet).
Starting from a fresh Sage tarball::

    ./sage -i openssl
    make ssl

And if you've already built Sage::

    ./sage -i openssl
    ./sage -f python2
    make ssl

The third line will rebuild all parts of Sage that depend on Python;
this can take a while.

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
directory to build Sage and its `HTML <http://en.wikipedia.org/wiki/HTML>`_
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
  run ``make doc-clean``, because the presence, or lack, of pictures
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

- ``make doc-clean`` removes several directories which are produced
  when building the documentation.

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

  for clean upstream tarballs, and it searches the directories

  - ``SAGE_SERVER/spkg/standard/``,
  - ``SAGE_SERVER/spkg/optional/``,
  - ``SAGE_SERVER/spkg/experimental/``,
  - ``SAGE_SERVER/spkg/archive/``

  for old-style Sage packages.


Here are some of the more commonly used variables affecting the build process:

- :envvar:`MAKE` - one useful setting for this variable when building Sage is
  ``MAKE='make -jNUM'`` to tell the ``make`` program to run ``NUM`` jobs in
  parallel when building.
  Note that not all Sage packages (e.g. ATLAS) support this variable.

  Some people advise using more jobs than there are CPU cores, at least if the
  system is not heavily loaded and has plenty of RAM; for example, a good
  setting for ``NUM`` might be between 1 and 1.5 times the number of cores.
  In addition, the ``-l`` option sets a load limit: ``MAKE='make -j4 -l5.5``,
  for example, tells ``make`` to try to use four jobs, but to not start more
  than one job if the system load average is above 5.5.
  See the manual page for GNU ``make``: `Command-line options
  <http://www.gnu.org/software/make/manual/make.html#Options-Summary>`_
  and `Parallel building
  <http://www.gnu.org/software/make/manual/make.html#Parallel>`_.

  .. warning::

      Some users on single-core OS X machines have reported problems when
      building Sage with ``MAKE='make -jNUM'`` with ``NUM`` greater than one.

- :envvar:`SAGE_NUM_THREADS` - if set to a number, then when building the
  documentation, parallel doctesting, or running ``sage -b``, use this many
  threads.
  If this is not set, then determine the number of threads using the value of
  the :envvar:`MAKE` (see above) or :envvar:`MAKEFLAGS` environment variables.
  If none of these specifies a number of jobs, use one thread (except for
  parallel testing: there we use a default of the number of CPU cores, with a
  maximum of 8 and a minimum of 2).

- :envvar:`SAGE_CHECK` - if set to ``yes``, then during the build process,
  or when installing packages manually,
  run the test suite for each package which has one.
  See also :envvar:`SAGE_CHECK_PACKAGES`.

- :envvar:`SAGE_CHECK_PACKAGES` - if :envvar:`SAGE_CHECK` is set to ``yes``,
  then the default behavior is to run test suites for all spkgs which contain
  them.
  If :envvar:`SAGE_CHECK_PACKAGES` is set, it should be a comma-separated list
  of strings of the form ``package-name`` or ``!package-name``.
  An entry ``package-name`` means to run the test suite for the named package
  regardless of the setting of :envvar:`SAGE_CHECK`.
  An entry ``!package-name`` means to skip its test suite.
  So if this is set to ``mpir,!python``, then always run the test suite for
  MPIR, but always skip the test suite for Python.

  .. note::

     As of this writing (April 2013, Sage 5.8), the test suite for the Python
     spkg fails on most platforms.
     So when this variable is empty or unset, Sage uses a default of
     ``!python``.

- :envvar:`SAGE64` - if set to ``yes``, then build a 64-bit binary on platforms
  which default to 32-bit, even though they can build 64-bit binaries.
  It adds the compiler flag ``-m64`` when compiling programs.
  The :envvar:`SAGE64` variable is mainly of use on OS X (pre 10.6), Solaris
  and OpenSolaris, though it will add the ``-m64`` flag on any operating
  system.
  If you are running Linux or version 10.6 or later of OS X on a 64-bit
  machine, then Sage will automatically build a 64-bit binary, so this
  variable does not need to be set.

- :envvar:`CFLAG64` - default value ``-m64``.
  If Sage detects that it should build a 64-bit binary, then it uses this flag
  when compiling C code.
  Modify it if necessary for your system and C compiler.
  This should not be necessary on most systems -- this flag will typically be
  set automatically, based on the setting of :envvar:`SAGE64`, for example.

- :envvar:`SAGE_INSTALL_GCC` - by default, Sage will automatically detect
  whether to install the `GNU Compiler Collection (GCC) <http://gcc.gnu.org/>`_
  package or not (depending on whether C, C++, and Fortran compilers are present
  and the versions of those compilers).
  Setting ``SAGE_INSTALL_GCC=yes`` will force Sage to install GCC.
  Setting ``SAGE_INSTALL_GCC=no`` will prevent Sage from installing GCC.

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

- :envvar:`SAGE_DEBUG` - controls debugging support.
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

- :envvar:`SAGE_PROFILE` - controls profiling support. If this is set
  to ``yes``, profiling support is enabled where possible. Note that
  Python-level profiling is always avaliable; This option enables
  profiling in Cython modules.

- :envvar:`SAGE_SPKG_INSTALL_DOCS` - if set to ``yes``, then install
  package-specific documentation to
  :file:`$SAGE_ROOT/local/share/doc/PACKAGE_NAME/` when an spkg is
  installed.
  This option may not be supported by all spkgs.
  Some spkgs might also assume that certain programs are available on the
  system (for example, ``latex`` or ``pdflatex``).

- :envvar:`SAGE_DOC_MATHJAX` - by default, any LaTeX code in Sage's
  documentation is processed by MathJax.
  If this variable is set to ``no``, then MathJax is not used -- instead,
  math is processed using LaTeX and converted by dvipng to image files,
  and then those files are included into the documentation.
  Typically, building the documentation using LaTeX and dvipng takes longer
  and uses more memory and disk space than using MathJax.

- :envvar:`SAGE_DOCBUILD_OPTS` - the value of this variable is passed as an
  argument to ``sage --docbuild all html`` or ``sage --docbuild all pdf`` when
  you run ``make``, ``make doc``, or ``make doc-pdf``.
  For example, you can add ``--no-plot`` to this variable to avoid building
  the graphics coming from the ``.. PLOT`` directive within the documentation,
  or you can add ``--no-tests`` to omit all "TESTS" blocks in the
  reference manual. Run ``sage --docbuild help`` to see the full list
  of options.

- :envvar:`SAGE_BUILD_DIR` - the default behavior is to build each spkg in a
  subdirectory of :file:`$SAGE_ROOT/local/var/tmp/sage/build/`; for
  example, build version 3.8.3.p12 of
  :file:`atlas` in the directory
  :file:`$SAGE_ROOT/local/var/tmp/sage/build/atlas-3.8.3.p12/`.
  If this variable is set, then build in
  :file:`$SAGE_BUILD_DIR/atlas-3.8.3.p12/` instead.
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

- :envvar:`SAGE_FAT_BINARY` - to build binaries that will run on the
  widest range of target CPUs set this variable to ``yes`` before
  building Sage. This does not make the binaries relocatable, it only
  avoids newer CPU instruction set extensions. For relocatable (=can
  be moved to a different directory) binaries, you must use
  https://github.com/sagemath/binary-pkg

Variables to set if you're trying to build Sage with an unusual setup, e.g.,
an unsupported machine or an unusual compiler:

- :envvar:`SAGE_PORT` - if you try to build Sage on a platform which is
  recognized as being unsupported (e.g. AIX, or HP-UX), or with a compiler
  which is unsupported (anything except GCC), you will see a message saying
  something like::

      You are attempting to build Sage on IBM's AIX operating system,
      which is not a supported platform for Sage yet. Things may or
      may not work. If you would like to help port Sage to AIX,
      please join the sage-devel discussion list -- see
      http://groups.google.com/group/sage-devel
      The Sage community would also appreciate any patches you submit.

      To get past this message and try building Sage anyway,
      export the variable SAGE_PORT to something non-empty.

  If this is case and you want to try to build Sage anyway, follow the
  directions: set :envvar:`SAGE_PORT` to something non-empty (and expect to
  run into problems).

Environment variables dealing with specific Sage packages:

- :envvar:`SAGE_MP_LIBRARY` - to use an alternative library in place of ``MPIR``
  for multiprecision integer arithmetic. Supported values are

    ``MPIR`` (default choice), ``GMP``.

- :envvar:`SAGE_ATLAS_ARCH` - if you are compiling ATLAS (in particular,
  if :envvar:`SAGE_ATLAS_LIB` is not set), you can use this environment
  variable to set a particular architecture and instruction set extension,
  to control the maximum number of threads ATLAS can use, and to trigger the
  installation of a static library (which is disabled by default unless
  building our custom shared libraries fails).
  The syntax is

    ``SAGE_ATLAS_ARCH=[threads:n,][static,]arch[,isaext1][,isaext2]...[,isaextN]``.

  While ATLAS comes with precomputed timings for a variety of CPUs, it only
  uses them if it finds an exact match.
  Otherwise, ATLAS runs through a lengthy automated tuning process in order
  to optimize performance for your particular system, which can take several
  days on slow and unusual systems.
  You drastically reduce the total Sage compile time if you manually select a
  suitable architecture.
  It is recommended to specify a suitable architecture on laptops or other
  systems with CPU throttling or if you want to distribute the binaries.
  Available architectures are

    ``POWER3``, ``POWER4``, ``POWER5``, ``PPCG4``, ``PPCG5``,
    ``POWER6``, ``POWER7``, ``IBMz9``, ``IBMz10``, ``IBMz196``,
    ``x86x87``, ``x86SSE1``, ``x86SSE2``, ``x86SSE3``, ``P5``,
    ``P5MMX``, ``PPRO``, ``PII``, ``PIII``, ``PM``, ``CoreSolo``,
    ``CoreDuo``, ``Core2Solo``, ``Core2``, ``Corei1``, ``Corei2``,
    ``Atom``, ``P4``, ``P4E``, ``Efficeon``, ``K7``, ``HAMMER``,
    ``AMD64K10h``, ``AMDDOZER``, ``UNKNOWNx86``, ``IA64Itan``,
    ``IA64Itan2``, ``USI``, ``USII``, ``USIII``, ``USIV``, ``UST2``,
    ``UnknownUS``, ``MIPSR1xK``, ``MIPSICE9``, ``ARMv7``.

  and instruction set extensions are

    ``VSX``, ``AltiVec``, ``AVXMAC``, ``AVXFMA4``, ``AVX``, ``SSE3``,
    ``SSE2``, ``SSE1``, ``3DNow``, ``NEON``.

  In addition, you can also set

  - ``SAGE_ATLAS_ARCH=fast`` which picks defaults for a modern (2-3 year old)
    CPU of your processor line, and

  - ``SAGE_ATLAS_ARCH=base`` which picks defaults that should work for a ~10
    year old CPU.

  For example,

    ``SAGE_ATLAS_ARCH=Corei2,AVX,SSE3,SSE2,SSE1``

  would be appropriate for a Core i7 CPU.

- :envvar:`SAGE_ATLAS_LIB` - if you have an installation of ATLAS on your
  system and you want Sage to use it instead of building and installing its
  own version of ATLAS, set this variable to be the directory containing your
  ATLAS installation.
  It should contain the files :file:`libatlas`, :file:`liblapack`,
  :file:`libcblas`, :file:`libf77blas` (and optionally :file:`libptcblas` and
  :file:`libptf77blas` for multi-threaded computations), with extensions ``.a``,
  ``.so``, or ``.dylib``.  For backward compatibility, the libraries may also be
  in the subdirectory :file:`SAGE_ATLAS_LIB/lib/`.

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

Sage uses the following environment variables when it runs:

- :envvar:`DOT_SAGE` - this is the directory, to which the user has read and
  write access, where Sage stores a number of files.
  The default location is :file:`$HOME/.sage/`.

- :envvar:`SAGE_STARTUP_FILE` - a file including commands to be executed every
  time Sage starts.
  The default value is :file:`$DOT_SAGE/init.sage`.

- :envvar:`SAGE_PATH` - a colon-separated list of directories which Sage
  searches when trying to locate Python libraries.

- :envvar:`SAGE_BROWSER` - on most platforms, Sage will detect the command to
  run a web browser, but if this doesn't seem to work on your machine, set this
  variable to the appropriate command.

- :envvar:`SAGE_CBLAS` - used in the file
  :file:`SAGE_ROOT/src/sage/misc/cython.py`.
  Set this to the base name of the BLAS library file on your system if you want
  to override the default setting.
  That is, if the relevant file is called :file:`libcblas_new.so` or
  :file:`libcblas_new.dylib`, then set this to ``cblas_new``.

Sage overrides the user's settings of the following variables:

- :envvar:`MPLCONFIGDIR` - ordinarily, this variable lets the user set their
  matplotlib config directory.
  Due to incompatibilies in the contents of this directory among different
  versions of matplotlib, Sage overrides the user's setting, defining it
  instead to be :file:`$DOT_SAGE/matplotlib-VER`, with ``VER`` replaced by the
  current matplotlib version number.

Variables dealing with doctesting:

- :envvar:`SAGE_TIMEOUT` - used for Sage's doctesting: the number of seconds
  to allow a doctest before timing it out.
  If this isn't set, the default is 300 seconds (5 minutes).

- :envvar:`SAGE_TIMEOUT_LONG` - used for Sage's doctesting: the number of
  seconds to allow a doctest before timing it out, if tests are run using
  ``sage -t --long``.
  If this isn't set, the default is 1800 seconds (30 minutes).

- :envvar:`SAGE_PICKLE_JAR` - if you want to update the the standard pickle
  jar, set this to something non-empty and run the doctest suite.
  See the documentation for the functions :func:`picklejar` and
  :func:`unpickle_all` in
  :file:`$SAGE_ROOT/src/sage/structure/sage_object.pyx`, online
  `here (picklejar)
  <http://doc.sagemath.org/html/en/reference/sage/structure/sage_object.html#sage.structure.sage_object.picklejar>`_
  and `here (unpickle_all)
  <http://doc.sagemath.org/html/en/reference/sage/structure/sage_object.html#sage.structure.sage_object.unpickle_all>`_.

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

System-wide install
~~~~~~~~~~~~~~~~~~~

In the instructions below, we assume that ``/path/to/sage-x.y`` is
the directory where you want to install Sage.

#. First of all, extract the Sage source tarball in ``/path/to``
   (this will create the directory ``/path/to/sage-x.y``).
   After extracting, you can change the directory name if you do not
   like ``sage-x.y``.

#. Change the ownership of the ``/path/to/sage-x.y`` directory tree
   to your normal user account (as opposed to ``root``). This is because
   Sage will refuse to compile as ``root``. ::

       chown -R user:group /path/to/sage-x.y

#. Using your normal user account, build Sage.
   See the :ref:`build-from-source-step-by-step` above.

#. Make a symbolic link to the ``sage`` script in :file:`/usr/local/bin`::

       ln -s /path/to/sage-x.y/sage /usr/local/bin/sage

   Alternatively, copy the Sage script::

       cp /path/to/sage-x.y/sage /usr/local/bin/sage

   If you do this, make sure you edit the line::

       #SAGE_ROOT=/path/to/sage-version

   at the beginning of the copied ``sage`` script according to the direction
   given there to something like::

       SAGE_ROOT=<SAGE_ROOT>

   (note that you have to change ``<SAGE_ROOT>`` above!).
   It is recommended not to edit the original ``sage`` script, only the copy at
   :file:`/usr/local/bin/sage`.

#. Optionally, you can test Sage by running::

       make testlong

   or ``make ptestlong`` which tests files in parallel using multiple
   processes.
   You can also omit ``long`` to skip tests which take a long time.

Some common problems
--------------------

ATLAS
~~~~~

Usually Sage will build ATLAS with architectural defaults that are not tuned
to your particular CPU.
In particular, if your CPU has powersaving enabled then no accurate timings
can be made to tune the ATLAS build for your hardware.
If BLAS performance is critical for you, you must recompile ATLAS after
installing Sage either with architecture settings that match your hardware,
or run through ATLAS' automatic tuning process where timings of different
implementations are compared and the best choice used to build a custom ATLAS
library.
To do so, you have to

- Leave the computer idle while you are reinstalling ATLAS.
  Most of ATLAS will intentionally only compile/run on a single core.
  Accurate timings of cache edges require that the CPU is otherwise idle.

- Make sure that CPU powersaving mode (that is, anything but the
  ``performance`` CPU scaling governor in Linux) is turned off when building
  ATLAS.
  This requires administrator privileges.

- If your architecture is listed in :envvar:`SAGE_ATLAS_ARCH`, you should set
  it as it can help ATLAS in narrowing down the timing search.

To help you disable CPU power saving, Sage includes an ``atlas-config`` script
that will turn off CPU powersave and rebuild ATLAS.
The script will call ``sudo`` to gain the necessary rights, which may prompt
you for your password. For example::

    atlas-config

will run through the full automated tuning, and::

    SAGE_ATLAS_ARCH=Corei2,AVX,SSE3,SSE2,SSE1 atlas-config

would be appropriate if you have a Core i3/5/7 processor with AVX support.



**This page was last updated in July 2016 (Sage 7.3).**
