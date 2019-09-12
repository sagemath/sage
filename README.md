<a href="https://sagemath.org"><img src="src/doc/common/themes/sage/static/logo_sagemath_black.svg" height="60" align="right" /></a>

#  Sage: Open Source Mathematical Software

>   "Creating a Viable Open Source Alternative to
>    Magma, Maple, Mathematica, and MATLAB"

>   Copyright (C) 2005-2018 The Sage Development Team

   https://www.sagemath.org

The Sage Library is GPLv2+, and included packages have [compatible OSS
licenses](./COPYING.txt). [Over 400 people](https://www.sagemath.org/development-map.html)
have contributed code to Sage. In many cases, documentation for modules
and functions list the authors.

Getting Started
---------------

If you downloaded a [binary](http://www.sagemath.org/download.html)
(i.e. a version of SageMath prepared for a specific operating system),
Sage is ready to start -- just open a terminal in the directory where
you extracted the binary archive and type:

    ./sage

If you downloaded the [sources](http://www.sagemath.org/download-source.html),
please read below on how to build Sage and work around common issues.

If you have questions or encounter problems, please do not hesitate
to email the [sage-support mailing list](https://groups.google.com/group/sage-support)
or ask on [ask.sagemath.org](https://ask.sagemath.org).

Contributing to Sage
--------------------

If you'd like to contribute to Sage, be sure to read the
[Developer's Guide](https://doc.sagemath.org/html/en/developer/index.html).

Supported Platforms
-------------------

Sage fully supports several Linux distributions, recent versions of
Mac OS X, Windows (using virtualization), as well as a number of
Solaris and OpenSolaris releases.

Ports are in progress to some other, less common platforms. The list of
supported platforms and their current statuses are given in [our wiki](https://wiki.sagemath.org/SupportedPlatforms).

If you are interested in helping port Sage to a new platform, please let
us know at the [sage-devel mailing list](https://groups.google.com/group/sage-devel).

Quick Instructions to Build from Source
---------------------------------------

The following steps briefly outline the process of building Sage from
source. More detailed instructions, including how to build faster on
multicore machines, are contained later in this README and in the
[Installation Guide](https://doc.sagemath.org/html/en/installation).

1. Make sure your system has an SSL library and its development
files installed

   Like Python, on which it is based, Sage uses the OpenSSL library
   for added performance if made available by the operating system. It
   has been shown that Sage can be successfully built against other
   SSL libraries, with some of its features disabled.


1. Make sure you have the dependencies and 5 GB of free disk space

   * __All Linux versions:__ gcc, make, m4, perl, ranlib, git, and tar (a
   matching set of gcc, gfortran and g++ will avoid the compilation
   of Sage-specific compilers). It should also be possible to use clang/clang++,
   however this is less well-tested.

   * __Fedora or RedHat systems:__ the perl-ExtUtils-MakeMaker package.
   (install these using your package manager)

   * __OS X:__
       * Make sure you have installed the most recent version
       of Xcode which you can install for free from the App Store.
       * You also need to install the "command line tools". When
       using OS X Mavericks, after installing Xcode, run
       `xcode-select --install` from a terminal window:
       Then click "Install" in the pop-up window.
       When using OS X Mountain Lion or earlier, you need to install the
       command line tools from Xcode: run Xcode; then from the File
       menu, choose "Preferences", then the "Downloads" tab, and then
       "Install" the Command Line Tools. You might also have Homebrew or
       a similar "Apple's missing package manager" system installed, with
       and libraries such gfortran, gmp, etc installed. (However, this
       is still experimental as of May 2019).

   * __Other platforms:__ See detailed instructions below.

1. It might be desirable, it terms of faster building and better portability,
   to install, as system packages, an ever increasing [list of Sage packages](https://trac.sagemath.org/ticket/27330)
   which otherwise might have to be built. The following is a list of Sage packages
   "replaceable" by system's packages as of Sage release 8.8:
   `bzip2`, `curl`, `cmake`, `gcc/clang`, `gf2x`, `gfortran` (usually part of `gcc` installation),
   `git`, `gmp`, `libffi`, `patch`, `pcre`, `perl_term_readline_gnu`, `xz/lzma`, `yasm`, `zeromq`,  `zlib`.
   Details and names of system packages containing these are system-dependent.  E.g. on Debian
   `bzip2` lives in `libbz2-dev`. More details on this are in Installation manual.

1. Extract the tarball

       tar zxvf sage-*.tar.gz

1. cd into the Sage directory and type make

       cd sage-*/
       make

   That's it! Everything is automatic and non-interactive. The build
   should work fine on all fully supported platforms. If it does not, we
   want to know!

Environment Variables
---------------------

There are a lot of environment variables which control the install
process of Sage described in more detail in the
[Installation Guide](https://doc.sagemath.org/html/en/installation/source.html#environment-variables).

Implementation
--------------

Sage has significant components written in the following languages:
C/C++, Python, Cython, Lisp, Fortran, and a bit of Perl. Lisp (ECL), Python, and Cython
are built as part of Sage.

Docker Images
-------------

You can also have a look at our Docker images to run Sage.
To use these images [install Docker](https://www.docker.com/community-edition#/download)
and follow the instructions on [our Docker Hub page](https://hub.docker.com/r/sagemath/sagemath/).

More Detailed Instructions to Build from Source
-----------------------------------------------

1. Make sure you have about 5 GB of free disk space.

1. Install build dependencies

   * __Linux:__ See quick instructions above.

   * __OS X:__ (a.k.a __MacOS__) Make sure you have a recent Xcode version.
   If you don't, go to https://developer.apple.com/,
   sign up, and download the free Xcode package. Usually, Xcode's command line
   tools suffice to build Sage, although several times new releases of Xcode broke this.
   Only OS X >= 10.4 is supported, and (as of May 2019) we only test Sage on OS X >= 10.6.

   * __Solaris and OpenSolaris:__ Building Sage on these platforms is more
   tricky than on Linux or OS X. For details on how to build Sage on
   these platforms, see [our wiki](https://wiki.sagemath.org/solaris) (outdated as of May 2019).

   * __Windows:__ [Download and install VirtualBox](https://www.virtualbox.org/wiki/Downloads),
   and then download the [Sage virtual appliance](https://wiki.sagemath.org/SageAppliance).

   * __NOTE:__ On some operating systems, it might be necessary to install
   gas/as, gld/ld, gnm/nm. On most platforms, these are automatically
   installed when you install the programs listed above.

1. Extract the Sage source tarball into a directory, making sure
   there are no spaces in the path to the resulting directory.

   Note that moving the directory after Sage has been built will
   require to build Sage again.

1. Change to the Sage directory using `cd`.

1. Optional: set some environment variables to customize the build.

   For example, the `MAKE` environment variable controls whether to run
   several jobs in parallel, while the `SAGE_CHECK` environment variable
   controls whether to perform more tests during the installation.  For
   an in-depth discussion of environment variables for building Sage, see
   [the installation guide](http://doc.sagemath.org/html/en/installation/source.html#environment-variables).

   On a machine with 4 processors, say, typing `export MAKE="make -j4"`
   will configure the build script to perform a parallel compilation of
   Sage using 4 jobs. You might even consider `-j5` or `-j6`, as
   building with more jobs than CPU cores can speed things up further.
   You might in addition pass a `-l` [load flag](https://www.gnu.org/software/make/manual/make.html#Options-Summary)
   to `make`: this sets a load limit, so for example if you execute
   `export MAKE="make -j4 -l5.5"` then "make" won't start more than one
   job at a time if the system load average is above 5.5, see
   the [make documentation](https://www.gnu.org/software/make/manual/make.html#Parallel).

   If you want to run the test suite for each individual Sage package
   as it gets installed, type `export SAGE_CHECK="yes"`. This will run
   each test suite, raising an error if any failure occurs. Python's
   test suite has been disabled by default, because it causes failures
   on most systems. To enable the Python test suite, set the environment
   variable `SAGE_CHECK_PACKAGES` to `python`.

1. To start the build, type `make`.

   Note: to build a Python3-based Sage, instead of typing `make`, type

       make configure
       ./configure --with-python=3
       make

   This will build Sage based on Python 3 rather than based on Python 2,
   which is still the default at this point. The resulting Sage mostly
   works well, though some features (less of them at each release!) are
   not yet ready for Python 3. The progress on this is tracked at
   [Sage Trac ticket 15530: Metaticket: Add support for python 3.6+](https://trac.sagemath.org/ticket/15530).

1. Wait about 20 minutes to 14 days, depending on your computer (it took
   about 2 weeks to build Sage on the T-Mobile G1 Android cell phone).

1. Type `./sage` to try it out.

1. Optional: Type `make ptestlong` to test all examples in the documentation
   (over 200,000 lines of input!) -- this takes from 10 minutes to
   several hours. Don't get too disturbed if there are 2 to 3 failures,
   but always feel free to email the section of `logs/ptestlong.log` that
   contains errors to the [sage-support mailing list](https://groups.google.com/group/sage-support).
   If there are numerous failures, there was a serious problem with your build.

   Note: if you built for Python 3, you can instead run `make ptest-python3`.

1. The HTML version of the [documentation](http://doc.sagemath.org/html/en/index.html)
   is built during the compilation process of Sage and resides in the directory
   `local/share/doc/sage/html/`.

1. Optional: If you want to build the PDF version of the documentation,
    run `make doc-pdf` (this requires LaTeX to be installed).

1. Optional: You might install optional packages of interest to you: type
   `./sage --optional` to get a list.

1. Optional: It is recommended that you have both LaTeX and the
   ImageMagick tools (e.g. the "convert" command) installed since some
   plotting functionality benefits from it.

1. Optional: Read this if you are intending to run a Sage notebook
   server for multiple users. For security (i.e., to run
   `notebook(secure=True)`) you want to access the server using the
   HTTPS protocol. First, install OpenSSL and the OpenSSL development
   headers on your system if they are not already installed. Then
   install pyOpenSSL by building Sage and then typing
   `./sage -i pyopenssl`.
   Note that this command requires internet access. Alternatively,
   `make ssl` builds Sage and installs pyOpenSSL.

Troubleshooting
---------------

If you have problems building Sage, check the Sage Installation Guide,
and also note the following. Each separate component of Sage is
contained in an spkg; these are stored in `build/pkgs/`. As each one
is built, a build log is stored in `logs/pkgs/`, so you can browse these
to find error messages. If an spkg fails to build, the whole build
process will stop soon after, so check the most recent log files
first, or run

       grep -li "^Error" logs/pkgs/*

from the top-level Sage directory to find log files with error
messages in them.  Send (a small part of) the relevant log file to the
[sage-devel mailing list](https://groups.google.com/group/sage-devel),
making sure to include at least some of the error messages; probably
someone there will have some helpful suggestions.

Supported Compilers
-------------------

Sage includes a GCC (_GNU Compiler Collection_) package. However,
it almost always better to use C, C++ and Fortran compilers
already available one the system. To force using  specific compilers,
set environment variables `CC`, `CXX`, and `FC` (for C, C++, and Fortran compilers,
respectively) to the desired values,
and run `./configure`. E.g. `CC=clang CXX=clang++ FC=gfortran ./configure`
will configure Sage to be built with Clang C/C++ compilers and Fortran
compiler gfortran.

It is determined automatically whether Sage's GCC package, or just its part containing
Fortran compiler `gfortran` needs to be installed. This can be overwritten
by running `./configure` with option `--without-system-gcc`.

There are some known problems with old assemblers, in particular when
building the ECM package. You should ensure that your assembler
understands all instructions for your processor. On Linux, this means
you need a recent version of binutils; on OS X you need a recent version
of Xcode.

Directory Layout
----------------

Simplified directory layout (only essential files/directories):
```
SAGE_ROOT                 Root directory (sage-x.y.z in Sage tarball)
├── build
│   ├── deps              Dependency information of packages
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
│   ├── *.m4
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
building and installing each of Sage's dependencies (e.g. `make python2`).

Although it's possible to manually run Sage's `configure` script if one wants
to provide some customizations (e.g. it is possible to select which BLAS
implementation to use), the top-level `Makefile` will run `configure` for you,
in order to build `build/make/Makefile` since it's a prerequisite for most of
Sage's make targets.

The `configure` script itself, if it is not already built, can be generated by
running the `bootstrap` script (the latter requires _GNU autotools_ being installed).
The top-level `Makefile` also takes care of this automatically.

To summarize, running a command like `make python2` at the top-level of the
source tree goes something like this:

1. `make python2`
1. run `./bootstrap` if `configure` does not exist
1. run `./configure` if `build/make/Makefile` does not exist
1. `cd` into `build/make` and run the `install` script--this is little more
   than a front-end to running `make -f build/make/Makefile python2`, which
   sets some necessary environment variables and logs some information
1. `build/make/Makefile` contains the actual rule for building `python2`; this
   includes building all of `python2`'s dependencies first (and their
   dependencies, recursively); the actual package installation is performed
   with the `sage-spkg` program


Relocation
----------

It used to be possible to move the `sage-x.y.z/` directory anywhere you
want, however, this is no longer supported.
If you copy the sage script or make a symbolic link to it, you
should modify the script to reflect this (as instructed at the top of
the script). It is important that the path to Sage does not have any spaces
and non-ASCII characters in it.

For a system-wide installation, you have to build Sage as a "normal" user
and then as root you can change permissions. Afterwards, you need to start up
Sage as root at least once prior to using the system-wide Sage as a
normal user. See the [Installation Guide](https://doc.sagemath.org/html/en/installation/source.html#installation-in-a-multiuser-environment)
for further information.

If you find anything that doesn't work correctly after you moved the
directory, please email the [sage-support mailing list](https://groups.google.com/group/sage-support).

Redistribution
--------------

Your local Sage install is almost exactly the same as any "developer"
install. You can make changes to documentation, source, etc., and very
easily package the complete results up for redistribution just like we
do.

1. To make your own source tarball of Sage, type:

       sage --sdist

   The result is placed in the directory `dist/`.

2. To make a binary distribution with your currently installed packages,
   visit [sagemath/binary-pkg](https://github.com/sagemath/binary-pkg).


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
