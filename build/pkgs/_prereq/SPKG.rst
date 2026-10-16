_prereq: Represents system packages required for installing SageMath from source
================================================================================

Description
-----------

This dummy package represents the minimal requirements (system packages)
for installing SageMath from source.

In addition to standard :wikipedia:`POSIX <POSIX>` utilities
and the :wikipedia:`bash <Bash_(Unix_shell)>` shell,
the following standard command-line development tools must be installed on your
computer:

- **C and C++ compilers**: GCC (GNU Compiler Collection) or Clang (LLVM).
  C99 and C++11 support are required. Use matching compiler vendors and versions.
  On macOS, use Apple's Clang compilers from the Xcode Command Line Tools.
- **Fortran compiler**: a compiler accepting free-form source code.
  Officially, we support GNU ``gfortran``; LLVM ``flang`` may also work.
  Supported GCC and gfortran versions are 10.3 through 16.x.
- **make**: GNU make, version 3.80 or later. Version 3.82 or later is recommended.
- **m4**: GNU m4 1.4.2 or later (non-GNU or older versions might also work).
- **perl**: version 5.8.0 or later.
- **ar** and **ranlib**: can be obtained as part of GNU binutils.
- **tar**: GNU tar version 1.17 or later, or BSD tar (as provided on macOS).
- **patch**.
- **boost**: the library ``boost`` with its headers.
- **bzip2**: the executable ``bzip2`` and the library ``libbz2`` with its headers.
  (some Linux distros package these separately, e.g. Debian/Ubuntu needs
  packages ``bzip2`` and ``libbz2-dev``; Fedora needs ``bzip2`` and ``bzip2-devel``.)
- **GNU readline/history**: the libraries and development headers, including
  ``readline/readline.h``. Apple's libedit compatibility library does not
  provide the required GNU API; install GNU readline separately on macOS.
- **ncurses/terminfo**: the terminal libraries and development headers, including
  ``ncurses.h``. A wide-character ncurses installation is also supported.
  Some distributions split the terminfo library into a separate package;
  install the dependencies of their readline and ncurses development packages.
- **pkgconf**, also known as ``pkg-config``.
- **zlib**: the library ``libz`` with its headers, and its pkg-config zlib.pc file.
  (some Linux distros package these separately.)
  On macOS we provide a basic zlib.pc file in build/platform/macos/pkgconfig,
  which, if needed, gets prepended to PKG_CONFIG_PATH by ./configure.

Other versions of these may work, but they are untested.

On macOS, suitable versions of most of these tools are provided
by the Xcode Command Line Tools.  To install them, after installing XCode itself,
open a terminal window and run ``xcode-select --install``; then click "Install" in the
pop-up window.  If the Xcode Command Line Tools are already installed,
you may want to check if they need to be updated by typing
``softwareupdate -l``. The remaining are provided by either one of macOS's
"missing package managers", such as Homebrew, or as standalone
tools. In particular ``pkgconf`` is available as `pkg-config_pkg
<https://github.com/donmccaughey/pkg-config_pkg>`_.

On Linux, ``ar`` and ``ranlib`` are in the `binutils
<https://www.gnu.org/software/binutils/>`_ package.  The other
programs are usually located in packages with their respective names.

Boost is available with most supported distribitions.
It can also be installed using the Boost's project `installer B2
<https://www.boost.org/doc/user-guide/getting-started.html#_download_boost>`_.
After downloading and untarring the archive, and changing to the directory with the sources,
``./bootstrap.sh && ./b2 && ./b2 install --prefix=/usr/local`` will
install Boost in ``/usr/local``; this takes around 5 minutes of wall clock time
on a moderately fast M1 Apple Mac. Instead of ``/usr/local`` one may choose another location,
say ``/opt/foo``,which then might have to be passed (in case the location is not known
to the toolchain) to Sage via its ``./configure``,
with ``--with-boost=/opt/foo`` option.

On Redhat-derived systems not all perl components are installed by
default and you might have to install the ``perl-ExtUtils-MakeMaker``
package.

To check if you have the above prerequisites installed, for example ``perl``,
type::

    $ command -v perl

or::

    $ which perl

on the command line. If it gives an error (or returns nothing), then
either ``perl`` is not installed, or it is installed but not in your
:wikipedia:`PATH <PATH_%28variable%29>`.

Selecting compilers
-------------------

Sage uses system compilers and does not build or install GCC or gfortran.
Select the C, C++, and Fortran compilers with the environment variables
:envvar:`CC`, :envvar:`CXX`, and :envvar:`FC` when running ``./configure``.
For example::

    $ ./configure CC=clang CXX=clang++ FC=gfortran

Users of older Linux distributions should upgrade before building Sage.
In particular, use Ubuntu 22.04 (Jammy) or newer. To select a matching
versioned GNU toolchain, for example::

    $ sudo apt-get install gcc-15 g++-15 gfortran-15
    $ ./configure CC=gcc-15 CXX=g++-15 FC=gfortran-15

If these packages are unavailable in your release's repositories, install
suitable compilers using your distribution's toolchain repositories or ask
your system administrator. On macOS, Homebrew's ``gcc`` package provides
``gfortran``; continue to use Apple's Clang for C and C++.
