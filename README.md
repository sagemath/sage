<a href="https://sagemath.org"><img src="src/doc/common/themes/sage/static/logo_sagemath_black.svg" height="60" align="right" /></a>

#  Sage: Open Source Mathematical Software

>   "Creating a Viable Open Source Alternative to
>    Magma, Maple, Mathematica, and MATLAB"

>   Copyright (C) 2005-2020 The Sage Development Team

   https://www.sagemath.org

The Sage Library is GPLv2+, and included packages have [compatible OSS
licenses](./COPYING.txt). [Over 400 people](https://www.sagemath.org/development-map.html)
have contributed code to Sage. In many cases, documentation for modules
and functions list the authors.

Getting Started
---------------

If you downloaded a [binary](https://www.sagemath.org/download.html)
(i.e. a version of SageMath prepared for a specific operating system),
Sage is ready to start -- just open a terminal in the directory where
you extracted the binary archive and type:

    ./sage

(Note that the first run will take more time, as Sage needs to get itself ready.)

If you downloaded the [sources](https://www.sagemath.org/download-source.html),
please read below on how to build Sage and work around common issues.

If you have questions or encounter problems, please do not hesitate
to email the [sage-support mailing list](https://groups.google.com/group/sage-support)
or ask on [ask.sagemath.org](https://ask.sagemath.org).

Supported Platforms
-------------------

Sage fully supports all major Linux distributions, recent versions of
macOS, and Windows (using Cygwin, Windows Subsystem for Linux, or
using virtualization).

We highly appreciate contributions to Sage that fix portability bugs
and help port Sage to new platforms; let us know at the [sage-devel
mailing list](https://groups.google.com/group/sage-devel).

Docker Images
-------------

You can also have a look at our Docker images to run Sage.
To use these images [install Docker](https://www.docker.com/community-edition#/download)
and follow the instructions on [our Docker Hub page](https://hub.docker.com/r/sagemath/sagemath/).

[Windows] Preparing the Platform
--------------------------------

The 64-bit version of Cygwin, also known as Cygwin64, is the current
target for Sage support on Windows.

1. Download [cygwin64](https://cygwin.com/install.html) (do not get
   the 32-bit version; it is not supported by Sage).

2. Run the `setup-x86_64.exe` graphical installer.  Pick the default
   options in most cases.  At the package selection screen, use the
   search bar to find and select at least the following packages:
   `bzip2`, `coreutils`, `curl`, `gawk`, `gzip`, `tar`, `wget`, `git`.

3. Start the Cygwin terminal and ensure you get a working bash prompt.

4. Make sure the path of your Cygwin home directory does not contain
   space characters.

   By default, your username in Cygwin is the same as your username in
   Windows.  This might contain spaces and other traditionally
   non-UNIX-friendly characters, e.g., if it is your full name.  You
   can check this as follows:

        $ whoami
        Erik M. Bray

   This means your default home directory on Cygwin contains this
   username verbatim; in the above example, `/home/Erik M. Bray`.  It
   will save some potential trouble if you change your Cygwin home
   directory to something not containing any non-alphanumeric
   characters, for example, `/home/embray`.  The easiest way to do
   this is to first create the home directory you want to use instead,
   then create an `/etc/passwd` file specifying that directory as your
   home, as follows:

        $ whocanibe=embray
        $ mkdir /home/$whocanibe
        $ mkpasswd.exe -l -u "$(whoami)" | sed -r 's,/home/[^:]+,/home/'$whocanibe, > /etc/passwd

   After this, close all Cygwin terminals (ensure nothing in
   `C:\cygwin64` is running), then start a new Cygwin terminal and
   your home directory should have moved.

   There are [other ways to do
   this](https://stackoverflow.com/questions/1494658/how-can-i-change-my-cygwin-home-folder-after-installation),
   but the above seems to be the simplest that's still supported.

5. Install the package manager `apt-cyg`:

        $ curl -OL https://rawgit.com/transcode-open/apt-cyg/master/apt-cyg
        $ install apt-cyg /usr/local/bin
        $ rm -f apt-cyg

An alternative to Cygwin is to use [Windows Subsystem for
Linux](https://docs.microsoft.com/en-us/windows/wsl/faq), which allows
you to install a standard Linux distribution such as Ubuntu within
your Windows.  Then all instructions for installation in Linux apply.

As another alternative, you can also run Linux on Windows using Docker
(see above) or other virtualization solutions such as the [Sage
virtual appliance](https://wiki.sagemath.org/SageAppliance).

[macOS] Preparing the Platform
------------------------------

Make sure you have installed the most current version of Xcode
supported on your version of macOS.  If you don't, go to
https://developer.apple.com/, sign up, and download the free Xcode
package.

You also need to install the "command line tools": After installing
Xcode, run `xcode-select --install` from a terminal window; then click
"Install" in the pop-up window.  (When using Mountain Lion or earlier,
you need to install the command line tools from Xcode: run Xcode; then
from the File menu, choose "Preferences", then the "Downloads" tab,
and then "Install" the Command Line Tools.)

Optionally, you can consider installing Homebrew ("the missing package
manager for macOS") from https://brew.sh/, which can provide libraries
such gfortran, gmp, etc.

Instructions to Build from Source
---------------------------------

Like many other software packages, Sage is built from source using
`./configure`, followed by `make`.  However, we strongly recommend to
read the following step-by-step instructions for building Sage.

The instructions cover all of Linux, macOS, and Cygwin.

More detailed instructions are contained in the [Installation
Guide](https://doc.sagemath.org/html/en/installation).

1. Decide on the source/build directory (`SAGE_ROOT`):

    - For example, you could use `SAGE_ROOT=~/sage/sage-x.y`, which we
      will use as the running example below, where `x.y` is the
      current Sage version.

    - You need at least 6 GB of free disk space.

    - The path name must contain **no spaces**.

    - After starting the build, you cannot move the source/build
      directory without breaking things.

    - [Cygwin] Avoid building in home directories of Windows domain
      users or in paths with capital letters.

2. Download/unpack the sources.

    - After downloading the source tarball `sage-x.y.tar.gz` into
      `~/sage/`:

            $ cd ~/sage/
            $ tar zxvf sage-x.y.tar.gz

      This will create the subdirectory `sage-x.y`.

    - [Git] Alternatively, clone the Sage git repository:

            $ git clone -c core.symlinks=true --branch master https://github.com/sagemath/sage.git

      This will create the subdirectory `sage`. `cd sage/` and pick the branch you need
      by doing `git checkout` - typically you want the latest development branch, thus do
      `git checkout develop`.

    - [Windows] The Sage source tree contains symbolic links, and the
      build will not work if Windows line endings rather than UNIX
      line endings are used.

      Therefore it is crucial that you unpack the source tree from the
      Cygwin (or WSL) `bash` using the Cygwin (or WSL) `tar` utility
      and not using other Windows tools (including mingw). Likewise,
      when using `git`, it is recommended (but not necessary) to use the Cygwin (or WSL)
      version of `git`.

3. `cd` into the source/build directory:

        $ cd sage*/

4. Optionally, decide on the installation prefix (`SAGE_LOCAL`):

    - Traditionally, and by default, Sage is installed into the
      subdirectory hierarchy rooted at `SAGE_ROOT/local`.

    - This can be changed using `./configure --prefix=SAGE_LOCAL`,
      where `SAGE_LOCAL` is the desired installation prefix, which
      must be writable by the user.  (See the installation manual for
      options if you want to install into shared locations such as
      `/usr/local/`.  Do not attempt to build Sage as root.)

5. [Git] If you cloned the Sage repository using `git`, bootstrap the
   source tree using:

        $ make configure

6. [Linux, Cygwin] Install the required minimal build prerequisites.

   * Compilers: `gcc`, `gfortran`, `g++` (a matching set of these
     three will avoid the compilation of Sage-specific compilers -
     unless they are too old).  See the Installation Manual for a
     discussion of suitable compilers.

   * Build tools: GNU `make`, GNU `m4`, `perl` (including
     ``ExtUtils::MakeMaker``), `ranlib`, `git`, `tar`, `bc`

   * Python 3.4 or later, or Python 2.7, a full installation including
     `urllib`; but ideally version 3.7.x, 3.8.x, or 3.9.x, which will
     avoid having to build Sage's own copy of Python 3.

   We have collected lists of system packages that provide these build
   prerequisites. See [build/pkgs/arch.txt](build/pkgs/arch.txt),
   [cygwin.txt](build/pkgs/cygwin.txt),
   [debian.txt](build/pkgs/debian.txt) (also for Ubuntu, Linux Mint,
   etc.), [fedora.txt](build/pkgs/fedora.txt) (also for Red Hat,
   CentOS), and [slackware.txt](build/pkgs/slackware.txt).

7. Optional, but highly recommended: Make sure your system has an SSL
   library and its development files installed.

   Like Python, on which it is based, Sage uses the OpenSSL library
   for added performance if made available by the operating system. It
   has been shown that Sage can be successfully built against other
   SSL libraries, with some of its features disabled.

8. Optional: It is recommended that you have both LaTeX and the
   ImageMagick tools (e.g. the "convert" command) installed since some
   plotting functionality benefits from it.

9. Optionally, review the configuration options, which includes
   many optional packages:

        ./configure --help

10. Optional, but highly recommended: Set some environment variables to
   customize the build.

   For example, the `MAKE` environment variable controls whether to
   run several jobs in parallel.  On a machine with 4 processors, say,
   typing `export MAKE="make -j4"` will configure the build script to
   perform a parallel compilation of Sage using 4 jobs. On some
   powerful machines, you might even consider `-j16`, as building with
   more jobs than CPU cores can speed things up further.

   To reduce the terminal output during the build, type `export V=0`.
   (`V` stands for "verbosity".)

   For an in-depth discussion of more environment variables for
   building Sage, see [the installation
   guide](https://doc.sagemath.org/html/en/installation/source.html#environment-variables).

11. Type `./configure`, followed by any options that you wish to use.
   For example, to build Sage with `gf2x` package supplied by Sage,
   use `./configure --with-system-gf2x=no`.

   At the end of a successful `./configure` run, you may see messages
   recommending to install extra system packages using your package
   manager.

   For a large [list of Sage
   packages](https://trac.sagemath.org/ticket/27330), Sage is able to
   detect whether an installed system package is suitable for use with
   Sage; in that case, Sage will not build another copy from source.

   Sometimes, the messages will recommend to install packages that are
   already installed on your system. See the earlier configure
   messages or the file `config.log` for explanation.  Also, the
   messages may recommend to install packages that are actually not
   available; only the most recent releases of your distribution will
   have all of these recommended packages.

12. Optional: If you choose to install the additional system packages,
   a re-run of `./configure` will test whether the versions installed
   are usable for Sage; if they are, this will reduce the compilation
   time and disk space needed by Sage. The usage of packages may be
   adjusted by `./configure` parameters (check again the output of
   `./configure --help`).

13. Type `make`.  That's it! Everything is automatic and
   non-interactive; but it will take a few hours (on a recent
   computer).

   The build should work fine on all fully supported platforms. If it
   does not, we want to know!

14. Type `./sage` to try it out.

15. Optional: Type `make ptestlong` to test all examples in the documentation
   (over 200,000 lines of input!) -- this takes from 10 minutes to
   several hours. Don't get too disturbed if there are 2 to 3 failures,
   but always feel free to email the section of `logs/ptestlong.log` that
   contains errors to the [sage-support mailing list](https://groups.google.com/group/sage-support).
   If there are numerous failures, there was a serious problem with your build.

16. The HTML version of the [documentation](https://doc.sagemath.org/html/en/index.html)
   is built during the compilation process of Sage and resides in the directory
   `local/share/doc/sage/html/`.

17. Optional: If you want to build the PDF version of the documentation,
    run `make doc-pdf` (this requires LaTeX to be installed).

18. Optional: You might install optional packages of interest to you: type
   `./sage --optional` to get a list.

Troubleshooting
---------------

If you have problems building Sage, check the Sage Installation Guide,
as well as the version-specific Sage Installation FAQ in the [Sage Release
Tour](https://wiki.sagemath.org/ReleaseTours) corresponding to the
version that you are installing.

Please do not hesitate to ask for help in the [SageMath forum
](https://ask.sagemath.org/questions/) or the [sage-support mailing
list](https://groups.google.com/forum/#!forum/sage-support).  The
[Troubleshooting section in the Sage Installation Guide
](https://doc.sagemath.org/html/en/installation/troubles.html)
provides instructions on what information to provide so that we can provide
help more effectively.

Contributing to Sage
--------------------

If you'd like to contribute to Sage, we strongly recommend that you read the
[Developer's Guide](https://doc.sagemath.org/html/en/developer/index.html).

Sage has significant components written in the following languages:
C/C++, Python, Cython, Common Lisp, Fortran, and a bit of Perl.

Directory Layout
----------------

Simplified directory layout (only essential files/directories):
```
SAGE_ROOT                 Root directory (sage-x.y.z in Sage tarball)
├── build
│   └── pkgs              Every package is a subdirectory here
│       ├── atlas
│       …
│       └── zn_poly
├── COPYING.txt           Copyright information
├── local                 Compiled packages are installed here
│   ├── bin               Executables
│   ├── include           C/C++ headers
│   ├── lib               Shared libraries
│   ├── share             Databases, architecture-independent data, docs
│       └── doc           Viewable docs of Sage and of some components
│   └── var
│       ├── sage          List of installed packages
│       └── tmp           Temporary files when building Sage
├── logs
│   ├── dochtml.log       Log of the documentation build
│   ├── install.log       Full install log
│   └── pkgs              Build logs of individual packages
│       ├── atlas-3.10.1.p7.log
│       …
│       └── zn_poly-0.9.p11.log
├── m4                    M4 macros for configure
│   └── *.m4
├── Makefile              Running "make" uses this file
├── README.md             This file
├── sage                  Script to start Sage
├── src                   All of Sage source (not third-party packages)
│   ├── bin               Scripts that Sage uses internally
│   ├── doc               Sage documentation sources
│   └── sage              The Sage library source code
├── upstream              Source tarballs of packages
│   ├── atlas-3.10.1.tar.bz2
│   …
│   └── zn_poly-0.9.tar.bz2
└── VERSION.txt
```
For more details see [our Developer's Guide](https://doc.sagemath.org/html/en/developer/coding_basics.html#files-and-directory-structure).


Build System
------------

This is a brief summary of the Sage software distribution's build system.
There are two components to the full Sage system--the Sage Python library
and its associated user interfaces, and the larger software distribution of
Sage's main dependencies (for those dependencies not supplied by the user's
system).

Sage's Python library is built and installed using a `setup.py` script as is
standard for Python packages (Sage's `setup.py` is non-trivial, but not
unusual).

Most of the rest of the build system is concerned with building all of Sage's
dependencies in the correct order in relation to each other.  The dependencies
included by Sage are referred to as SPKGs (i.e. "Sage Packages") and are listed
under `build/pkgs`.

The main entrypoint to Sage's build system is the top-level `Makefile` at the
root of the source tree.  Unlike most normal projects that use autoconf (Sage
does as well, as described below), this `Makefile` is not generated.  Instead,
it contains a few high-level targets and targets related to bootstrapping the
system.  Nonetheless, we still run `make <target>` from the root of the source
tree--targets not explicitly defined in the top-level `Makefile` are passed
through to another Makefile under `build/make/Makefile`.

The latter `build/make/Makefile` *is* generated by an autoconf-generated
`configure` script, using the template in `build/make/Makefile.in`.  This
includes rules for building the Sage library itself (`make sagelib`), and for
building and installing each of Sage's dependencies (e.g. `make gf2x`).

The `configure` script itself, if it is not already built, can be generated by
running the `bootstrap` script (the latter requires _GNU autotools_ being installed).
The top-level `Makefile` also takes care of this automatically.

To summarize, running a command like `make python3` at the top-level of the
source tree goes something like this:

1. `make python3`
2. run `./bootstrap` if `configure` needs updating
3. run `./configure` with any previously configured options if `build/make/Makefile`
   needs updating
4. `cd` into `build/make` and run the `install` script--this is little more
   than a front-end to running `make -f build/make/Makefile python3`, which
   sets some necessary environment variables and logs some information
5. `build/make/Makefile` contains the actual rule for building `python3`; this
   includes building all of `python3`'s dependencies first (and their
   dependencies, recursively); the actual package installation is performed
   with the `sage-spkg` program


Relocation
----------

It is not supported to move the `SAGE_ROOT` or `SAGE_LOCAL` directory
after starting the build.  If you do move the directories, you will
have to rebuilt Sage again from scratch.

If you copy the `sage` script or make a symbolic link to it, you
should modify the script to reflect this (as instructed at the top of
the script). It is important that the path to Sage does not have any
spaces and non-ASCII characters in it.

For a system-wide installation, you have to build Sage as a "normal" user
and then as root you can change permissions. Afterwards, you need to start up
Sage as root at least once prior to using the system-wide Sage as a
normal user. See the [Installation Guide](https://doc.sagemath.org/html/en/installation/source.html#installation-in-a-multiuser-environment)
for further information.

Redistribution
--------------

Your local Sage install is almost exactly the same as any "developer"
install. You can make changes to documentation, source, etc., and very
easily package the complete results up for redistribution just like we
do.

1. To make a binary distribution with your currently installed packages,
   visit [sagemath/binary-pkg](https://github.com/sagemath/binary-pkg).

2. (**Obsolete, probably broken**) To make your own source tarball of Sage, type:

       sage --sdist

   The result is placed in the directory `dist/`.


Changes to Included Software
----------------------------

All software included with Sage is copyrighted by the respective authors
and released under an open source license that is __GPL version 3 or
later__ compatible. See [COPYING.txt](./COPYING.txt) for more details.

Sources are in unmodified (as far as possible) tarballs in the
`upstream/` directory. The remaining description, version
information, patches, and build scripts are in the accompanying
`build/pkgs/<packagename>` directory. This directory is
part of the Sage git repository.
