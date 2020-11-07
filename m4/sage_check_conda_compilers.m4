AC_DEFUN([SAGE_CHECK_CONDA_COMPILERS], [
    AC_MSG_CHECKING([whether a conda environment is active])
    AS_IF([test "x$CONDA_PREFIX" != x], [have_conda_active=yes], [have_conda_active=no])
    AC_MSG_RESULT($have_conda_active)
    AS_IF([test $have_conda_active = yes], [
        dnl A conda environment is active.
        dnl #27699: Conda compiler packages must be installed
        dnl #30662: Do not allow system pkg-config to give us library paths that bypass conda
        AS_VAR_SET([reasons], [])
        AS_VAR_SET([need_pkgs], [])
        AS_IF([test -z "$CC" -o -z "$CXX" -o -z "$FC" ], [
            AS_VAR_APPEND([reasons],
                          [", but
at least one of the environment variables CC, CXX, FC is not set"])
            AS_VAR_APPEND([need_pkgs],
                          [" c-compiler cxx-compiler fortran-compiler"])
        ])
        AS_CASE(["$PKG_CONFIG"],
          [$CONDA_PREFIX/bin/*pkg-config], [
            dnl pkg-config from conda, possibly with target prefix
          ], [
            dnl system pkg-config, reject
            AS_IF([test -z "$reasons"],
                  [AS_VAR_APPEND([reasons], [", but"])],
                  [AS_VAR_APPEND([reasons], [" and"])])
            AS_VAR_APPEND([reasons],
                          ["
the pkg-config command is not provided by the conda environment"])
            AS_VAR_APPEND([need_pkgs],
                          [" pkg-config"])
        ])
        AS_IF([test -n "$need_pkgs"], [
            AC_MSG_ERROR([A conda environment ($CONDA_DEFAULT_ENV) is active$reasons,
which indicates that the conda environment is missing
the following conda packages required for building Sage:
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
