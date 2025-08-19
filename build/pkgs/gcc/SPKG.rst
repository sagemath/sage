gcc: The GNU Compiler Collection or other suitable C and C++ compilers
======================================================================

Description
-----------

This package represents the required C and C++ compilers.

- GCC (GNU Compiler Collection) versions 8.x (>= 8.4.0) to 13.x are supported.

- Clang (LLVM) is also supported.

The required Fortran compiler is represented by the package ``gfortran``.

You can pass the names of compilers to use to ``./configure`` using
the environment variables :envvar:`CC`, :envvar:`CXX`, and
:envvar:`FC`, for C, C++, and Fortran compilers, respectively.

For example, if your C compiler is ``clang``, your C++ compiler is
``clang++``, and your Fortran compiler is ``flang``, then you would
need to run::

    $ ./configure CC=clang CXX=clang++ FC=flang

Vendor and versions of the C and C++ compilers should match.

Users of older Linux distributions (in particular, ``ubuntu-xenial``
or older, ``debian-buster`` or older, ``linuxmint-18`` or older)
should upgrade their systems before attempting to install Sage from
source.  Users of ``ubuntu-bionic``, ``linuxmint-19.x``, and
``opensuse-15.x`` can install a versioned ``gcc`` system package
and then use::

    $ ./configure CC=gcc-8 CXX=g++-8 FC=gfortran-8

or similar. Users on ``ubuntu`` can also install a modern compiler
toolchain `using the ubuntu-toolchain-r ppa
<https://askubuntu.com/questions/1140183/install-gcc-9-on-ubuntu-18-04/1149383#1149383>`_.
On ``ubuntu-trusty``, also the package ``binutils-2.26`` is required;
after installing it, make it available using ``export
PATH="/usr/lib/binutils-2.26/bin:$PATH"``.  Instead of upgrading their
distribution, users of ``centos-7`` can install a modern compiler
toolchain `using Redhat's devtoolset
<https://stackoverflow.com/a/67212990/557937>`_.

This package uses the non-standard default
``configure --with-system-gcc=force``, giving an error at ``configure``
time when no suitable system compilers are configured.

You can override this using ``./configure --without-system-gcc``.  In
this case, Sage builds and installs the GNU Compiler Collection,
including the C, C++ and Fortran compiler. This is not recommended.
You will need suitable C and C++ compilers from which GCC can
bootstrap itself. There are some known problems with old assemblers,
in particular when building the ``ecm`` and ``fflas_ffpack``
packages. You should ensure that your assembler understands all
instructions for your processor. On Linux, this means you need a
recent version of ``binutils`` (not provided by an SPKG); on macOS
you need a recent version of Xcode.

(Installing the
``gfortran`` SPKG becomes a no-op in this case.)

Building Sage from source on Apple Silicon (M1, M2, M3, M4; arm64) requires
the use of Apple's Command Line Tools, and those tools include a suitable
compiler. Sage's ``gcc`` SPKG is not suitable for Apple Silicon; building it
will likely fail.

License
-------

GPL version 2 or version 3


Upstream Contact
----------------

https://gcc.gnu.org/
