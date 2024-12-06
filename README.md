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

Flag algebras
-----------------

This repository is a copy of the official SageMath project with additional functionality to handle 
flag algebraic calculations. 

Table of Contents
-----------------

* [Instructions to Build from Source](#instructions-to-build-from-source)
* [Troubleshooting](#troubleshooting)
* [Directory Layout](#directory-layout)
* [Build System](#build-system)
* [Relocation](#relocation)
* [Redistribution](#redistribution)
* [Changes to Included Software](#changes-to-included-software)


Instructions to Build from Source
---------------------------------

Sage attempts to support all major Linux distributions, recent versions of
macOS, and Windows (using Windows Subsystem for Linux or
virtualization). The additional software and packages needed for flag algebraic calculations
are only tested on Linux. This guide only describes the installation process for Linux.

More details, providing a background for these instructions, can be found
in the section [Install from Source Code](https://doc.sagemath.org/html/en/installation/source.html)
in the Installation Guide.


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

3.  Install system packages.

    Either refer for this to the [section on installation from
    source](https://doc.sagemath.org/html/en/installation/source.html) in the
    Sage Installation Manual for compilations of system packages
    that you can install. When done, skip to step 7 (bootstrapping).

    Alternatively, follow the more fine-grained approach below.

4.  [Linux, WSL] Install the required minimal build prerequisites:

    - Compilers: `gcc`, `gfortran`, `g++` (GCC versions from 8.4.0 to 13.x
      and recent versions of Clang (LLVM) are supported).
      See [build/pkgs/gcc/SPKG.rst](build/pkgs/gcc/SPKG.rst) and
      [build/pkgs/gfortran/SPKG.rst](build/pkgs/gfortran/SPKG.rst)
      for a discussion of suitable compilers.

    - Build tools: GNU `make`, GNU `m4`, `perl` (including
      `ExtUtils::MakeMaker`), `ranlib`, `git`, `tar`, `bc`.
      See [build/pkgs/_prereq/SPKG.rst](build/pkgs/_prereq/SPKG.rst) for
      more details.

    - Python 3.4 or later, or Python 2.7, a full installation including
      `urllib`; but ideally version 3.9.x, 3.10.x, 3.11.x, 3.12.x, which
      will avoid having to build Sage's own copy of Python 3.
      See [build/pkgs/python3/SPKG.rst](build/pkgs/python3/SPKG.rst)
      for more details.

    We have collected lists of system packages that provide these build
    prerequisites. See, in the folder
    [build/pkgs/_prereq/distros](build/pkgs/_prereq/distros),
    the files
    [arch.txt](build/pkgs/_prereq/distros/arch.txt),
    [debian.txt](build/pkgs/_prereq/distros/debian.txt)
    (also for Ubuntu, Linux Mint, etc.),
    [fedora.txt](build/pkgs/_prereq/distros/fedora.txt)
    (also for Red Hat, CentOS),
    [opensuse.txt](build/pkgs/_prereq/distros/opensuse.txt),
    [slackware.txt](build/pkgs/_prereq/distros/slackware.txt), and
    [void.txt](build/pkgs/_prereq/distros/void.txt), or visit
    https://doc.sagemath.org/html/en/reference/spkg/_prereq.html#spkg-prereq

5.  Optional: It is recommended that you have both LaTeX and
    the ImageMagick tools (e.g. the "convert" command) installed
    since some plotting functionality benefits from them.

6.  [Development] If you plan to do Sage development or otherwise work with
    ticket branches and not only releases, install the bootstrapping
    prerequisites. See the files in the folder
    [build/pkgs/_bootstrap/distros](build/pkgs/_bootstrap/distros), or
    visit
    https://doc.sagemath.org/html/en/reference/spkg/_bootstrap.html#spkg-bootstrap

7.  Bootstrap the source tree using the following command:

        $ make configure

    (If the bootstrapping prerequisites are not installed, this command
    will download a package providing pre-built bootstrap output instead.)

8.  Sanitize the build environment. Use the command

        $ env

    to inspect the current environment variables, in particular `PATH`,
    `PKG_CONFIG_PATH`, `LD_LIBRARY_PATH`, `CFLAGS`, `CPPFLAGS`, `CXXFLAGS`,
    and `LDFLAGS` (if set).

    Remove items from these (colon-separated) environment variables
    that Sage should not use for its own build. In particular, remove
    items if they refer to a previous Sage installation.

    - [WSL] In particular, WSL imports many items from the Windows
      `PATH` variable into the Linux environment, which can lead to
      confusing build errors. These items typically start with `/mnt/c`.
      It is best to remove all of them from the environment variables.
      For example, you can set `PATH` using the command:

            $ export PATH=/usr/sbin/:/sbin/:/bin/:/usr/lib/wsl/lib/

    - [macOS with homebrew] Set required environment variables for the build:

            $ source ./.homebrew-build-env

      This is to make some of Homebrew's packages (so-called keg-only
      packages) available for the build. Run it once to apply the
      suggestions for the current terminal session. You may need to
      repeat this command before you rebuild Sage from a new terminal
      session, or after installing additional homebrew packages.  (You
      can also add it to your shell profile so that it gets run
      automatically in all future sessions.)

9.  Optionally, decide on the installation prefix (`SAGE_LOCAL`):

    - Traditionally, and by default, Sage is installed into the
      subdirectory hierarchy rooted at `SAGE_ROOT/local/`.

    - This can be changed using `./configure --prefix=SAGE_LOCAL`,
      where `SAGE_LOCAL` is the desired installation prefix, which
      must be writable by the user.

      If you use this option in combination with `--disable-editable`,
      you can delete the entire Sage source tree after completing
      the build process.  What is installed in `SAGE_LOCAL` will be
      a self-contained installation of Sage.

    - Note that in Sage's build process, `make` builds **and**
      installs (`make install` is a no-op).  Therefore the
      installation hierarchy must be writable by the user.

    - See the Sage Installation Manual for options if you want to
      install into shared locations such as `/usr/local/`.
      Do not attempt to build Sage as `root`.

10. Optionally, review the configuration options, which includes
    many optional packages:

        $ ./configure --help

    Notable options for Sage developers are the following:

    - Use the option `--config-cache` to have `configure`
      keep a disk cache of configuration values. This gives a nice speedup
      when trying out ticket branches that make package upgrades, which
      involves automatic re-runs of the configuration step.

    - Use the option `--enable-ccache` to have Sage install and use the
      optional package `ccache`, which is preconfigured to keep a
      disk cache of object files created from source files. This can give
      a great speedup when switching between different branches, at the
      expense of disk space use.

11. Optional, but highly recommended: Set some environment variables to
    customize the build.

    For example, the `MAKE` environment variable controls whether to
    run several jobs in parallel.  On a machine with 4 processors, say,
    typing `export MAKE="make -j4"` will configure the build script to
    perform a parallel compilation of Sage using 4 jobs. On some
    powerful machines, you might even consider `-j16`, as building with
    more jobs than CPU cores can speed things up further.

    To reduce the terminal output during the build, type `export V=0`.
    (`V` stands for "verbosity".)

    Some environment variables deserve a special mention: `CC`,
    `CXX` and `FC`. These variables defining your compilers
    can be set at configuration time and their values will be recorded for
    further use at build time and runtime.

    For an in-depth discussion of more environment variables for
    building Sage, see [the installation
    guide](https://doc.sagemath.org/html/en/installation/source.html#environment-variables).

12. Type `./configure`, followed by any options that you wish to use.
    For example, to build Sage with `gf2x` package supplied by Sage,
    use `./configure --with-system-gf2x=no`.

    At the end of a successful `./configure` run, you may see messages
    recommending to install extra system packages using your package
    manager.

    For a large [list of Sage
    packages](https://github.com/sagemath/sage/issues/27330), Sage is able to
    detect whether an installed system package is suitable for use with
    Sage; in that case, Sage will not build another copy from source.

    Sometimes, the messages will recommend to install packages that are
    already installed on your system. See the earlier configure
    messages or the file `config.log` for explanation.  Also, the
    messages may recommend to install packages that are actually not
    available; only the most recent releases of your distribution will
    have all of these recommended packages.

13. Optional: If you choose to install the additional system packages,
    a re-run of `./configure` will test whether the versions installed
    are usable for Sage; if they are, this will reduce the compilation
    time and disk space needed by Sage. The usage of packages may be
    adjusted by `./configure` parameters (check again the output of
    `./configure --help`).

14. Type `make`.  That's it! Everything is automatic and
    non-interactive.

    If you followed the above instructions, in particular regarding the
    installation of system packages recommended by the output of
    `./configure` (step 11), and regarding the parallel build (step 10),
    building Sage takes less than one hour on a modern computer.
    (Otherwise, it can take much longer.)

    The build should work fine on all fully supported platforms. If it
    does not, we want to know!

15. Type `./sage` to try it out. In Sage, try for example `2 + 2`,
    `plot(x^2)`, `plot3d(lambda x, y: x*y, (-1, 1), (-1, 1))`
    to test a simple computation and plotting in 2D and 3D.
    Type <kbd>Ctrl</kbd>+<kbd>D</kbd> or `quit` to quit Sage.

16. Optional: Type `make ptestlong` to test all examples in the documentation
    (over 200,000 lines of input!) -- this takes from 10 minutes to
    several hours. Don't get too disturbed if there are 2 to 3 failures,
    but always feel free to email the section of `logs/ptestlong.log` that
    contains errors to the [sage-support mailing list](https://groups.google.com/group/sage-support).
    If there are numerous failures, there was a serious problem with your build.

17. The HTML version of the [documentation](https://doc.sagemath.org/html/en/index.html)
    is built during the compilation process of Sage and resides in the directory
    `local/share/doc/sage/html/`. You may want to bookmark it in your browser.

18. Optional: If you want to build the PDF version of the documentation,
    run `make doc-pdf` (this requires LaTeX to be installed).

19. Optional: Install optional packages of interest to you:
    get a list by typing  `./sage --optional` or by visiting the
    [packages documentation page](https://doc.sagemath.org/html/en/reference/spkg/).

20. Optional: Create a symlink to the installed `sage` script in a
    directory in your `PATH`, for example `/usr/local`. This will
    allow you to start Sage by typing `sage` from anywhere rather than
    having to either type the full path or navigate to the Sage
    directory and type `./sage`. This can be done by running:

        $ sudo ln -s $(./sage -sh -c 'ls $SAGE_ROOT/venv/bin/sage') /usr/local/bin

21. Optional: Set up SageMath as a Jupyter kernel in an existing Jupyter notebook
    or JupyterLab installation, as described in the section
    [Launching SageMath](https://doc.sagemath.org/html/en/installation/launching.html)
    in the Sage Installation Guide.

Alternative Installation using PyPI
---------------

For installing Sage in a Python environment from PyPI, Sage provides the
`pip`-installable package [sagemath-standard](https://pypi.org/project/sagemath-standard/).

Unless you need to install Sage into a specific existing environment, we recommend
to create and activate a fresh virtual environment, for example `~/sage-venv/`:

            $ python3 -m venv ~/sage-venv
            $ source ~/sage-venv/bin/activate

As the first installation step, install [sage_conf](https://pypi.org/project/sage-conf/),
which builds various prerequisite packages in a subdirectory of `~/.sage/`:

            (sage-venv) $ python3 -m pip install -v sage_conf

After a successful installation, a wheelhouse provides various Python packages.
You can list the wheels using the command:

            (sage-venv) $ ls $(sage-config SAGE_SPKG_WHEELS)

If this gives an error saying that `sage-config` is not found, check any messages
that the `pip install` command may have printed. You may need to adjust your `PATH`,
for example by:

            $ export PATH="$(python3 -c 'import sysconfig; print(sysconfig.get_path("scripts", "posix_user"))'):$PATH"

Now install the packages from the wheelhouse and the [sage_setup](https://pypi.org/project/sage-conf/)
package, and finally install the Sage library:

            (sage-venv) $ python3 -m pip install $(sage-config SAGE_SPKG_WHEELS)/*.whl sage_setup
            (sage-venv) $ python3 -m pip install --no-build-isolation -v sagemath-standard

The above instructions install the latest stable release of Sage.
To install the latest development version instead, add the switch `--pre` to all invocations of
`python3 -m pip install`.

**NOTE:** PyPI has various other `pip`-installable packages with the word "sage" in their names.
Some of them are maintained by the SageMath project, some are provided by SageMath users for
various purposes, and others are entirely unrelated to SageMath. Do not use the packages
`sage` and `sagemath`. For a curated list of packages, see the chapter
[Packages and Features](https://doc.sagemath.org/html/en/reference/spkg/index.html) of the
Sage Reference Manual.

Troubleshooting
---------------

If you have problems building Sage, check the Sage Installation Guide,
as well as the version-specific installation help in the [release
tour](https://github.com/sagemath/sage/releases) corresponding to the
version that you are installing.

Please do not hesitate to ask for help in the [SageMath forum
](https://ask.sagemath.org/questions/) or the [sage-support mailing
list](https://groups.google.com/forum/#!forum/sage-support).  The
[Troubleshooting section in the Sage Installation Guide
](https://doc.sagemath.org/html/en/installation/troubles.html)
provides instructions on what information to provide so that we can provide
help more effectively.

Directory Layout
----------------

Simplified directory layout (only essential files/directories):
```
SAGE_ROOT                 Root directory (create by git clone)
├── build
│   └── pkgs              Every package is a subdirectory here
│       ├── 4ti2/
│       …
│       └── zlib/
├── configure             Top-level configure script
├── COPYING.txt           Copyright information
├── pkgs                  Source trees of Python distribution packages
│   ├── sage-conf
│   │   ├── sage_conf.py
│   │   └── setup.py
│   ├── sage-docbuild
│   │   ├── sage_docbuild/
│   │   └── setup.py
│   ├── sage-setup
│   │   ├── sage_setup/
│   │   └── setup.py
│   ├── sage-sws2rst
│   │   ├── sage_sws2rst/
│   │   └── setup.py
│   └── sagemath-standard
│       ├── bin/
│       ├── sage -> ../../src/sage
│       └── setup.py
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
│       │   └── venv-python3.9  (SAGE_VENV)
│       │       │         Installation hierarchy (virtual environment)
│       │       │         for Python packages
│       │       ├── bin/  Executables and installed scripts
│       │       ├── lib/python3.9/site-packages/
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
│       └── zlib-1.2.11.log
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
│   └── zlib-1.2.11.tar.gz
├── venv -> SAGE_VENV     Convenience symlink to the virtual environment
└── VERSION.txt
```

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

Relocation
----------

It is not supported to move the `SAGE_ROOT` or `SAGE_LOCAL` directory
after building Sage.  If you do move the directories, you will have to
run ``make distclean`` and build Sage again from scratch.

For a system-wide installation, you have to build Sage as a "normal" user
and then as root you can change permissions. See the [Installation Guide](https://doc.sagemath.org/html/en/installation/source.html#installation-in-a-multiuser-environment)
for further information.

Redistribution
--------------

Your local Sage install is almost exactly the same as any "developer"
install. You can make changes to documentation, source, etc., and very
easily package the complete results up for redistribution just like we
do.

1.  To make a binary distribution with your currently installed packages,
    visit [sagemath/binary-pkg](https://github.com/sagemath/binary-pkg).

2.  To make your own source tarball of Sage, type:

        $ make dist

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

<p align="center">
   Copyright (C) 2005-2024 The Sage Development Team
</p>
<p align="center">
   https://www.sagemath.org
</p>

