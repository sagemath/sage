<div>
<a href="https://sagemath.org">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="src/doc/common/static/logo_sagemath_white.svg">
    <img src="src/doc/common/static/logo_sagemath_black.svg" height="60" align="left">
  </picture>
</a>
   <em>"Creating a Viable Open Source Alternative to
   Magma, Maple, Mathematica, and MATLAB"</em>
</div>

#

Sage is open source mathematical software released under the GNU General Public
Licence GPLv2+, and includes packages that have [compatible software licenses](./COPYING.txt).
[People all around the globe](https://www.sagemath.org/development-map.html) have contributed to the
development of Sage. [Full documentation](https://doc.sagemath.org/html/en/index.html) is available online.

## Table of Contents

- [Getting Started](#getting-started)
- [Supported Platforms](#supported-platforms)
- [\[Windows\] Preparing the Platform](#windows-preparing-the-platform)
- [\[macOS\] Preparing the Platform](#macos-preparing-the-platform)
- [Environment Setup](#environment-setup)
- [Instructions to Build from Source](#instructions-to-build-from-source)
- [SageMath Docker Images](#sagemath-docker-images)
- [Troubleshooting](#troubleshooting)
- [Contributing to Sage](#contributing-to-sage)
- [Directory Layout](#directory-layout)
- [Build System](#build-system)
- [Relocation](#relocation)
- [Redistribution](#redistribution)
- [Build System](#build-system)
- [Changes to Included Software](#changes-to-included-software)

## Getting Started

Those who are impatient may use prebuilt Sage available online from any of

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sagemath/sage-binder-env/master) &nbsp; [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-Ready--to--Code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/sagemath/sage/tree/master) &nbsp; [![Open in GitHub Codespaces](https://img.shields.io/badge/Open_in_GitHub_Codespaces-black?logo=github)](https://codespaces.new/sagemath/sage/tree/master)

without local installation. Otherwise read on.

The [Sage Installation Guide](https://doc.sagemath.org/html/en/installation/index.html)
provides a decision tree that guides you to the type of installation
that will work best for you. This includes building from source,
obtaining Sage from a package manager, using a container image, or using
Sage in the cloud.

If you have questions or encounter problems, please do not hesitate
to email the [sage-support mailing list](https://groups.google.com/group/sage-support)
or ask on the [Ask Sage questions and answers site](https://ask.sagemath.org).

## Supported Platforms

Sage attempts to support all major Linux distributions, recent versions of
macOS, and Windows (using Windows Subsystem for Linux or
virtualization).

Detailed information on supported platforms for a specific version of Sage
can be found in the section _Availability and installation help_ of the
[release tour for this version](https://github.com/sagemath/sage/releases).

We highly appreciate contributions to Sage that fix portability bugs
and help port Sage to new platforms; let us know at the [sage-devel
mailing list](https://groups.google.com/group/sage-devel).

## [Windows] Preparing the Platform

The preferred way to run Sage on Windows is using Windows Subsystem for
Linux (WSL). Follow the
[official WSL setup guide](https://docs.microsoft.com/en-us/windows/wsl/faq)
to install Ubuntu (or another Linux distribution).
Make sure you allocate WSL sufficient RAM; 5GB is known to work, while
2GB might be not enough for building Sage from source.
Then all instructions for installation in Linux apply.

As an alternative, you can also run Linux on Windows using Docker ([see
below](#sagemath-docker-images)) or other virtualization solutions.

## [macOS] Preparing the Platform

- If your Mac uses the Apple Silicon (M1, M2, M3, M4; arm64) architecture and
  you set up your Mac by transferring files from an older Mac, make sure
  that the directory `/usr/local` does not contain an old copy of Homebrew
  (or other software) for the x86_64 architecture that you may have copied
  over. Note that Homebrew for the M1 is installed in `/opt/homebrew`, not
  `/usr/local`.

- If you wish to use conda, please see the [section on
  conda](https://doc.sagemath.org/html/en/installation/conda.html) in the Sage
  Installation Manual for guidance.

- Otherwise, we strongly recommend to use Homebrew ("the missing package
  manager for macOS") from https://brew.sh/, which provides the `gfortran`
  compiler and many libraries.

- Otherwise, if you do not wish to install Homebrew, you will need to install
  the latest version of Xcode Command Line Tools. Open a terminal window and
  run `xcode-select --install`; then click "Install" in the pop-up window. If
  the Xcode Command Line Tools are already installed, you may want to check if
  they need to be updated by typing `softwareupdate -l`.

## Environment Setup

1.  Decide on the source/build directory (`SAGE_ROOT`):
    - On personal computers, any subdirectory of your :envvar:`HOME`
      directory should do.

    - For example, you could use `SAGE_ROOT=~/sage/sage`, which we
      will use as the running example below.

    - You need at least 10 GB of free disk space.

    - The full path to the source directory must contain **no spaces**.

    - After starting the build, you cannot move the source/build
      directory without breaking things.

    - You may want to avoid slow filesystems such as
      [network file systems (NFS)](https://en.wikipedia.org/wiki/Network_File_System)
      and the like.

    - [macOS] macOS allows changing directories without using exact capitalization.
      Beware of this convenience when compiling for macOS. Ignoring exact
      capitalization when changing into :envvar:`SAGE_ROOT` can lead to build
      errors for dependencies requiring exact capitalization in path names.

2.  Clone the sources with `git`:
    - To check that `git` is available, open a terminal and enter
      the following command at the shell prompt (`$`):

            $ git --version
            git version 2.42.0

      The exact version does not matter, but if this command gives an error,
      install `git` using your package manager, using one of these commands:

            $ sudo pacman -S git                          # on Arch Linux
            $ sudo apt-get update && apt-get install git  # on Debian/Ubuntu
            $ sudo yum install git                        # on Fedora/Redhat/CentOS
            $ sudo zypper install git                     # on openSUSE
            $ sudo xbps-install git                       # on Void Linux

    - Create the directory where `SAGE_ROOT` should be established:

            $ mkdir -p ~/sage
            $ cd ~/sage

    - Clone the Sage git repository:

            $ git clone -c core.symlinks=true --filter blob:none  \
                        --origin upstream --branch develop --tags \
                        https://github.com/sagemath/sage.git

      This command obtains the most recent development release.
      Replace `--branch develop` by `--branch master` to select
      the most recent stable release instead.

      This will create the subdirectory `~/sage/sage`. (See the section
      [Setting up git](https://doc.sagemath.org/html/en/developer/git_setup.html)
      and the following sections in the Sage Developer's Guide
      for more information.)

    - Change into the created subdirectory:

            $ cd sage

    - [Windows] The Sage source tree contains symbolic links, and the
      build will not work if Windows line endings rather than UNIX
      line endings are used.

      Therefore it is recommended (but not necessary) to use the
      WSL version of `git`.

## Instructions to Build from Source

Follow instructions in in the section [Install from Source Code](https://doc.sagemath.org/html/en/installation/source.html) in the Installation Guide.

## SageMath Docker images

[![Docker Status](http://dockeri.co/image/sagemath/sagemath)](https://hub.docker.com/r/sagemath/sagemath)

SageMath is available on Docker Hub and can be downloaded by:

```bash
docker pull sagemath/sagemath
```

Currently, only stable versions are kept up to date.

## Troubleshooting

If you have problems building Sage, check the Sage Installation Guide,
as well as the version-specific installation help in the [release
tour](https://github.com/sagemath/sage/releases) corresponding to the
version that you are installing.

Please do not hesitate to ask for help in the [SageMath forum
](https://ask.sagemath.org/questions/) or the [sage-support mailing
list](https://groups.google.com/forum/#!forum/sage-support). The
[Troubleshooting section in the Sage Installation Guide
](https://doc.sagemath.org/html/en/installation/troubles.html)
provides instructions on what information to provide so that we can provide
help more effectively.

## Contributing to Sage

If you'd like to contribute to Sage, we strongly recommend that you read the
[Developer's Guide](https://doc.sagemath.org/html/en/developer/index.html).

Sage has significant components written in the following languages:
C/C++, Python, Cython, Common Lisp, Fortran, and a bit of Perl.

## Directory Layout

Simplified directory layout (only essential files/directories):

```
SAGE_ROOT                 Root directory (create by git clone)
├── build
│   └── pkgs              Every package is a subdirectory here
│       ├── 4ti2/
│       …
│       └── zipp/
├── configure             Top-level configure script
├── COPYING.txt           Copyright information
├── local  (SAGE_LOCAL)   Installation hierarchy for non-Python packages
│   ├── bin               Executables
│   ├── include           C/C++ headers
│   ├── lib               Shared libraries, architecture-dependent data
│   ├── share             Databases, architecture-independent data, docs
│   │   └── doc           Viewable docs of Sage and of some components
│   └── var
│       ├── lib/sage
│       │   ├── installed/
│       │   │             Records of installed non-Python packages
│       │   ├── scripts/  Scripts for uninstalling installed packages
│       │   └── venv-python  (SAGE_VENV)
│       │       │         Installation hierarchy (virtual environment)
│       │       │         for Python packages
│       │       ├── bin/  Executables and installed scripts
│       │       ├── lib/python/site-packages/
│       │       │         Python modules/packages are installed here
│       │       └── var/lib/sage/
│       │           └── wheels/
│       │                 Python wheels for all installed Python packages
│       │
│       └── tmp/sage/     Temporary files when building Sage
├── logs
│   ├── install.log       Full install log
│   └── pkgs              Build logs of individual packages
│       ├── alabaster-0.7.12.log
│       …
│       └── zipp-3.19.0.log
├── m4                    M4 macros for generating the configure script
│   └── *.m4
├── Makefile              Running "make" uses this file
├── prefix -> SAGE_LOCAL  Convenience symlink to the installation tree
├── README.md             This file
├── sage                  Script to start Sage
├── src                   Monolithic Sage library source tree
│   ├── bin/              Scripts that Sage uses internally
│   ├── doc/              Sage documentation sources
│   └── sage/             The Sage library source code
├── upstream              Source tarballs of packages
│   ├── Babel-2.9.1.tar.gz
│   …
│   └── zipp-3.19.0-py3-none-any.whl
├── venv -> SAGE_VENV     Convenience symlink to the virtual environment
└── VERSION.txt
```

For more details see [our Developer's Guide](https://doc.sagemath.org/html/en/developer/coding_basics.html#files-and-directory-structure).

## Build System

This is a brief summary of the Sage software distribution's build system.
There are two components to the full Sage system--the Sage Python library
and its associated user interfaces, and the larger software distribution of
Sage's main dependencies (for those dependencies not supplied by the user's
system).

Sage's Python library is built and installed using a `setup.py` script as is
standard for Python packages (Sage's `setup.py` is non-trivial, but not
unusual).

Most of the rest of the build system is concerned with building all of Sage's
dependencies in the correct order in relation to each other. The dependencies
included by Sage are referred to as SPKGs (i.e. "Sage Packages") and are listed
under `build/pkgs`.

The main entrypoint to Sage's build system is the top-level `Makefile` at the
root of the source tree. Unlike most normal projects that use autoconf (Sage
does as well, as described below), this `Makefile` is not generated. Instead,
it contains a few high-level targets and targets related to bootstrapping the
system. Nonetheless, we still run `make <target>` from the root of the source
tree--targets not explicitly defined in the top-level `Makefile` are passed
through to another Makefile under `build/make/Makefile`.

The latter `build/make/Makefile` _is_ generated by an autoconf-generated
`configure` script, using the template in `build/make/Makefile.in`. This
includes rules for building the Sage library itself (`make sagelib`), and for
building and installing each of Sage's dependencies (e.g. `make gf2x`).

The `configure` script itself, if it is not already built, can be generated by
running the `bootstrap` script (the latter requires _GNU autotools_ being installed).
The top-level `Makefile` also takes care of this automatically.

To summarize, running a command like `make python3` at the top-level of the
source tree goes something like this:

1.  `make python3`
2.  run `./bootstrap` if `configure` needs updating
3.  run `./configure` with any previously configured options if `build/make/Makefile`
    needs updating
4.  change directory into `build/make` and run the `install` script--this is
    little more than a front-end to running `make -f build/make/Makefile python3`,
    which sets some necessary environment variables and logs some information
5.  `build/make/Makefile` contains the actual rule for building `python3`; this
    includes building all of `python3`'s dependencies first (and their
    dependencies, recursively); the actual package installation is performed
    with the `sage-spkg` program

## Relocation

It is not supported to move the `SAGE_ROOT` or `SAGE_LOCAL` directory
after building Sage. If you do move the directories, you will have to
run `make distclean` and build Sage again from scratch.

For a system-wide installation, you have to build Sage as a "normal" user
and then as root you can change permissions. See the [Installation Guide](https://doc.sagemath.org/html/en/installation/source.html#installation-in-a-multiuser-environment)
for further information.

## Redistribution

Your local Sage install is almost exactly the same as any "developer"
install. You can make changes to documentation, source, etc., and very
easily package the complete results up for redistribution just like we
do.

1.  To make a binary distribution with your currently installed packages,
    visit [sagemath/binary-pkg](https://github.com/sagemath/binary-pkg).

2.  To make your own source tarball of Sage, type:

        $ make dist

    The result is placed in the directory `dist/`.

## Changes to Included Software

All software included with Sage is copyrighted by the respective authors
and released under an open source license that is **GPL version 3 or
later** compatible. See [COPYING.txt](./COPYING.txt) for more details.

Sources are in unmodified (as far as possible) tarballs in the
`upstream/` directory. The remaining description, version
information, patches, and build scripts are in the accompanying
`build/pkgs/<packagename>` directory. This directory is
part of the Sage git repository.

<p align="center">
   Copyright (C) 2005-2025 The Sage Development Team
</p>
<p align="center">
   https://www.sagemath.org
</p>
