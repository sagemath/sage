gfortran: Fortran compiler from the GNU Compiler Collection
===========================================================

Description
-----------

This package represents the required Fortran compiler.

Officially we support ``gfortran`` from `GNU Compiler Collection (GCC)
<https://gcc.gnu.org/>`_.  It has also been reported that using ``flang``
(from LLVM) might work.

You can pass the names of compilers to use to ``./configure`` using
the environment variables :envvar:`CC`, :envvar:`CXX`, and
:envvar:`FC`, for C, C++, and Fortran compilers, respectively.

For example, if your C compiler is ``clang``, your C++ compiler is
``clang++``, and your Fortran compiler is ``flang``, then you would
need to run::

    $ ./configure CC=clang CXX=clang++ FC=flang

License
-------

GPL version 2 or version 3


Upstream Contact
----------------

http://gcc.gnu.org/

Special Update/Build Instructions
---------------------------------

None.
