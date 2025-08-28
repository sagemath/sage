SAGE_SPKG_CONFIGURE([openblas], [dnl CHECK
  SAGE_SPKG_DEPCHECK([gfortran], [dnl
    SAVE_LIBS="$LIBS"
    SAVE_CFLAGS="$CFLAGS"
    m4_pushdef([SAGE_OPENBLAS_MIN_VERSION_MAJOR], [0])
    m4_pushdef([SAGE_OPENBLAS_MIN_VERSION_MINOR], [2])
    m4_pushdef([SAGE_OPENBLAS_MIN_VERSION_MICRO], [20])
    m4_pushdef([SAGE_OPENBLAS_MIN_VERSION], [SAGE_OPENBLAS_MIN_VERSION_MAJOR.SAGE_OPENBLAS_MIN_VERSION_MINOR.SAGE_OPENBLAS_MIN_VERSION_MICRO])
    dnl Reject openblas 0.3.22 - https://github.com/sagemath/sage/pull/35371
    m4_pushdef([SAGE_OPENBLAS_LT_VERSION_MAJOR], [0])
    m4_pushdef([SAGE_OPENBLAS_LT_VERSION_MINOR], [3])
    m4_pushdef([SAGE_OPENBLAS_LT_VERSION_MICRO], [99])
    m4_pushdef([SAGE_OPENBLAS_LT_VERSION], [SAGE_OPENBLAS_LT_VERSION_MAJOR.SAGE_OPENBLAS_LT_VERSION_MINOR.SAGE_OPENBLAS_LT_VERSION_MICRO])
    PKG_CHECK_MODULES([OPENBLAS], [openblas >= ]SAGE_OPENBLAS_MIN_VERSION [openblas < ]SAGE_OPENBLAS_LT_VERSION, [dnl Have openblas.pc
      LIBS="$OPENBLAS_LIBS $LIBS"
      CFLAGS="$OPENBLAS_CFLAGS $CFLAGS"
      PKG_CHECK_VAR([OPENBLASPCDIR], [openblas], [pcfiledir], [dnl
        sage_install_blas_pc=yes
        AC_CHECK_FUNC([cblas_dgemm], [dnl openblas works as cblas
          sage_install_cblas_pc=yes
        ], [dnl openblas does not work as cblas; try to use system cblas as is
          PKG_CHECK_MODULES([CBLAS], [cblas], [], [sage_spkg_install_openblas=yes])
        ])
        dnl Check all name manglings that AC_FC_FUNC could check based on the
        dnl characteristics of the Fortran compiler
        m4_foreach([dgeqrf_mangled], [dgeqrf, dgeqrf_, DGEQRF, DGEQRF_], [dnl
          AC_CHECK_FUNC(dgeqrf_mangled, [dnl
            AS_VAR_SET([HAVE_DGEQRF], [yes])
          ])
        ])
        AS_IF([test x$HAVE_DGEQRF = xyes], [dnl openblas works as lapack
          sage_install_lapack_pc=yes
        ], [dnl openblas does not work as lapack; try to use system lapack as is
          PKG_CHECK_MODULES([LAPACK], [lapack], [], [sage_spkg_install_openblas=yes])
        ])
      ], [dnl
        AC_MSG_WARN([Unable to locate the directory of openblas.pc. This should not happen!])
        sage_spkg_install_openblas=yes
      ])
      AS_IF([test x$sage_spkg_install_openblas != xyes], [dnl
        AC_MSG_CHECKING([the OpenBLAS version using openblas_get_config])
        AC_LANG_PUSH([C])
        AC_RUN_IFELSE([dnl Reject 0.3.22 - see https://github.com/sagemath/sage/pull/35377
          AC_LANG_PROGRAM([[#include <string.h>
                            char *openblas_get_config(void); ]],
                          [[if (!strncmp(openblas_get_config(), "OpenBLAS 0.3.22", 15)) return 1;]])
        ], [dnl
          AC_MSG_RESULT([good])
        ], [dnl
          AC_MSG_RESULT([known bad version])
          sage_spkg_install_openblas=yes])
        AC_LANG_POP([C])
      ])
      AS_IF([test x$sage_spkg_install_openblas != xyes], [dnl
        AC_SUBST([SAGE_SYSTEM_FACADE_PC_FILES])
        AC_SUBST([SAGE_OPENBLAS_PC_COMMAND], ["\$(LN) -sf \"$OPENBLASPCDIR/openblas.pc\" \"\$(@)\""])
        m4_foreach([blaslibnam], [blas, cblas, lapack], [dnl
          AS_IF([test x$sage_install_]blaslibnam[_pc = xyes], [dnl
             AS_VAR_APPEND([SAGE_SYSTEM_FACADE_PC_FILES], [" \$(SAGE_PKGCONFIG)/]blaslibnam[.pc"])
          ])
        ])
      ])
    ], [dnl No openblas.pc
      dnl Recent OpenBLAS (>= 0.3.4, Dec 2018) provides the version number as
      dnl part of openblas_get_config.  We reject all older versions.
      AC_SEARCH_LIBS([openblas_get_config], [openblas cblas blas], [dnl
        AS_IF([test x"$ac_cv_search_openblas_get_config" != x"none required"], [dnl
          AS_VAR_APPEND([OPENBLAS_LIBS], ["$ac_cv_search_openblas_get_config "])
        ])
        AC_MSG_CHECKING([whether openblas_get_config indicates version >= ]SAGE_OPENBLAS_MIN_VERSION)
        AC_LANG_PUSH([C])
        AC_RUN_IFELSE([dnl
          AC_LANG_PROGRAM([[#include <stdio.h>
                            #include <string.h>
                            char *openblas_get_config(void);
                            int version[3]; ]],
                          [[version[0] = version[1] = version[2] = 0;
                            /*printf("%s", openblas_get_config());*/
                            if (sscanf(openblas_get_config(), "OpenBLAS %d.%d.%d", 
                                       version, version+1, version+2) < 1)
                              return 1;
                            if (  10000 * version[0]
                                  + 100 * version[1]
                                        + version[2]
                                < 10000 * ]]SAGE_OPENBLAS_MIN_VERSION_MAJOR[[
                                  + 100 * ]]SAGE_OPENBLAS_MIN_VERSION_MINOR[[
                                        + ]]SAGE_OPENBLAS_MIN_VERSION_MICRO[[)
                              return 1;
                            if (  10000 * version[0]
                                  + 100 * version[1]
                                        + version[2]
                                >=10000 * ]]SAGE_OPENBLAS_LT_VERSION_MAJOR[[
                                  + 100 * ]]SAGE_OPENBLAS_LT_VERSION_MINOR[[
                                        + ]]SAGE_OPENBLAS_LT_VERSION_MICRO[[)
                              return 1;
                            if (!strncmp(openblas_get_config(), "OpenBLAS 0.3.22", 15)) return 1;]])
        ], [AS_VAR_SET([HAVE_OPENBLAS], [yes])],
           [AS_VAR_SET([HAVE_OPENBLAS], [no])],
           [AS_VAR_SET([HAVE_OPENBLAS], [yes])])
        AC_LANG_POP([C])
        AC_MSG_RESULT([$HAVE_OPENBLAS])
      ])
      AC_SEARCH_LIBS([cblas_dgemm], [openblas cblas blas], [dnl
        AS_VAR_SET([HAVE_CBLAS_DGEMM], [yes])
        AS_IF([test x"$ac_cv_search_cblas_dgemm" != x"none required"], [dnl
          AS_VAR_APPEND([OPENBLAS_LIBS], ["$ac_cv_search_cblas_dgemm "])
        ])
      ], [], [-lgfortran])
      m4_foreach([dgeqrf_mangled], [dgeqrf, dgeqrf_, DGEQRF, DGEQRF_], [dnl
        AC_SEARCH_LIBS(dgeqrf_mangled, [openblas lapack], [dnl
          AS_VAR_SET([HAVE_DGEQRF], [yes])
          AS_IF([test x"$ac_cv_search_]dgeqrf_mangled[" != x"none required"], [dnl
            AS_VAR_APPEND([OPENBLAS_LIBS], ["$ac_cv_search_]dgeqrf_mangled[ "])
          ])
        ], [], [-lgfortran])
      ])
      AS_IF([test x"$HAVE_OPENBLAS" = xyes -a x"$HAVE_CBLAS_DGEMM" = xyes -a x"$HAVE_DGEQRF" = xyes], [dnl
        AC_SUBST([OPENBLAS_LIBS])
        AC_SUBST([SAGE_SYSTEM_FACADE_PC_FILES])
        AC_SUBST([SAGE_OPENBLAS_PC_COMMAND], ["  (echo \"Name: openblas\"; echo \"Description: OpenBLAS\"; echo \"Version: 0.3\"; echo \"Libs: $OPENBLAS_LIBS\") > \"\$(@)\""])
        m4_foreach([blaslibnam], [openblas, blas, cblas, lapack], [dnl
          AS_VAR_APPEND([SAGE_SYSTEM_FACADE_PC_FILES], [" \$(SAGE_PKGCONFIG)/]blaslibnam[.pc"])
        ])
      ], [dnl No system BLAS found
        sage_spkg_install_openblas=yes
      ])
    ])
    LIBS="$SAVE_LIBS"
    CFLAGS="$SAVE_CFLAGS"
  ])
])
