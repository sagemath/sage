AC_DEFUN([SAGE_CHECK_CONDA_COMPILERS], [
    AC_MSG_CHECKING([whether a conda environment is active])
    AS_IF([test "x$CONDA_PREFIX" != x], [have_conda_active=yes], [have_conda_active=no])
    AC_MSG_RESULT($have_conda_active)
    AS_IF([test $have_conda_active = yes], [
        dnl A conda environment is active.
        dnl #27699: Conda compiler packages must be installed
        need_pkgs="c-compiler cxx-compiler fortran-compiler"
        AS_IF([test -z "$CC" -o -z "$CXX" -o -z "$FC" ], [
          AC_MSG_ERROR([A conda environment ($CONDA_DEFAULT_ENV) is active, but
at least one of the environment variables CC, CXX, FC is not set, which indicates
that the conda environment is missing the following conda packages required
for building Sage:
    $need_pkgs
For building Sage, either:
- activate a conda environment that has these packages, using:
    conda activate ENVIRONMENT
- or install these conda packages, using
    conda install $need_pkgs
- or deactivate conda by
    conda deactivate
  (this command may need to be repeated).])
        ])
    ])
])
