gcc: The GNU Compiler Collection or other suitable C and C++ compilers
======================================================================

Description
-----------

This package represents the required C and C++ compilers.

- Sage's classical build checks for a C compiler that can compile C99 code.

- It also checks for a C++ compiler with C++11 support.

- GCC (GNU Compiler Collection) and Clang (LLVM) are both supported.

The required Fortran compiler is represented by the package ``gfortran``.

You can pass the names of compilers to use to ``./configure`` using
the environment variables :envvar:`CC`, :envvar:`CXX`, and
:envvar:`FC`, for C, C++, and Fortran compilers, respectively.

For example, if your C compiler is ``clang``, your C++ compiler is
``clang++``, and your Fortran compiler is ``flang``, then you would
need to run::

    $ ./configure CC=clang CXX=clang++ FC=flang

Vendor and versions of the C and C++ compilers should match.

Users of older Linux distributions should upgrade their systems before
attempting to install Sage from source.  In particular, users on
``ubuntu`` should use ``ubuntu-jammy`` (22.04) or a newer release.

The minimum supported GCC version is 10.3.  The following example uses
matching version ``15`` C, C++, and Fortran compilers.  On
``ubuntu-jammy``, these packages are available from
``ppa:ubuntu-toolchain-r/test``:

.. code-block:: bash

    $ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
    $ sudo apt-get update
    $ sudo apt-get install gcc-15 g++-15 gfortran-15

If these packages are already available from the standard repositories
of your release, omit the ``add-apt-repository`` command.  After
installation, select the compilers explicitly when configuring Sage::

    $ ./configure CC=gcc-15 CXX=g++-15 FC=gfortran-15

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
