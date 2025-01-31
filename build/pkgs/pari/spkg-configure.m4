SAGE_SPKG_CONFIGURE([pari], [
  dnl See gp_version below on how the version is computed from MAJV.MINV.PATCHV
  m4_pushdef([SAGE_PARI_MINVER],["134916"])dnl this version and higher allowed
  dnl Do not allow Pari 2.17 or later, see #38769:
  m4_pushdef([SAGE_PARI_MAXVER],["135424"])dnl this version and higher not allowed
  SAGE_SPKG_DEPCHECK([gmp readline], [
    AC_PATH_PROG([GP], [gp])
    if test x$GP = x; then dnl GP test
        AC_MSG_NOTICE([gp is not found])
        sage_spkg_install_pari=yes
    else
        AC_PATH_PROG([GPHELP], [gphelp])
        dnl needed for cypari2 installation; see #29319
        if test x$GPHELP = x; then
            AC_MSG_NOTICE([gphelp is not found; cannot use system pari/GP without gphelp])
            AC_MSG_NOTICE([Install a system package that provides it, possibly pari-doc.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari=yes
        else
            AC_MSG_CHECKING([whether gphelp has access to the documentation])
            dnl this is needed for cypari2, see #29342
            if $GPHELP -raw Catalan > /dev/null 2>&1; then
                AC_MSG_RESULT([yes])
            else
                AC_MSG_RESULT([no])
                AC_MSG_NOTICE([Install a system package that provides the documentation.])
                AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
                sage_spkg_install_pari=yes
            fi
        fi
        AC_MSG_CHECKING([is pari_elldata installed? ])
        gp_ell_check=`echo "r=ellinit(\"11a1\"); r[[11]]" | $GP -qf 2>> config.log`
        if test x$gp_ell_check = x20008; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without elldata package])
            AC_MSG_NOTICE([Install elldata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari=yes
        fi
        AC_MSG_CHECKING([is pari_galdata installed? ])
        gp_gal_check=`echo "polgalois(x^8-2)[[1]]" | $GP -qf 2>> config.log`
        if test x$gp_gal_check = x16; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without galdata package])
            AC_MSG_NOTICE([Install galdata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari=yes
        fi
        AC_MSG_CHECKING([is pari_galpol installed? ])
        gp_galp_check=`echo "galoisgetname(12,1) == \"C3 : C4\"" |  $GP -qf 2>> config.log`
        if test x$gp_galp_check = x1; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without galpol package])
            AC_MSG_NOTICE([Install galpol package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari=yes
        fi
        AC_MSG_CHECKING([is pari_seadata installed? ])
        gp_seadat_check=`echo "poldegree(ellmodulareqn(211)[[1]])" | $GP -qf 2>> config.log`
        if test x$gp_seadat_check = x212; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no; cannot use system pari/GP without seadata package])
            AC_MSG_NOTICE([Install seadata package and reconfigure.])
            AC_MSG_NOTICE([Otherwise Sage will build its own pari/GP.])
            sage_spkg_install_pari=yes
        fi

    fi dnl end GP test

      if test x$sage_spkg_install_pari = xno; then dnl main PARI test
        AC_CHECK_HEADER([pari/pari.h], [], [sage_spkg_install_pari=yes])
        dnl matpermanent appears in pari 2.11
        AC_SEARCH_LIBS([matpermanent], [pari], [
              AC_MSG_CHECKING([getting GP's version ])
              gp_version=`echo "v=version(); v[[1]]<<16 + v[[2]]<<8 + v[[3]]" | $GP -qf 2>> config.log`
              AC_MSG_RESULT([$gp_version])
              AC_MSG_CHECKING([comparing GP and libpari versions])
              AC_LANG_PUSH(C)
              AC_RUN_IFELSE([AC_LANG_PROGRAM([
                  [#include <pari/pari.h>]],
                 [[long vers;]
                  [pari_init(5000000, 2);]
                  [vers=paricfg_version_code;]
                  [pari_close()];
                  [return vers!=$gp_version;]])],
                [AC_MSG_RESULT([libpari's and GP's versions match. Good])],
                [AC_MSG_RESULT([libpari's version does not match GP's version. Not good])
                         sage_spkg_install_pari=yes],
                 [AC_MSG_RESULT([cross compiling. Assume they match])])
              AC_MSG_CHECKING([is GP's version good enough? ])
              AX_COMPARE_VERSION([$gp_version], [ge], [$SAGE_PARI_MINVER], [
               AX_COMPARE_VERSION([$gp_version], [lt], [$SAGE_PARI_MAXVER], [
                  AC_MSG_RESULT([yes])
                  AC_MSG_CHECKING([getting GP's datadir])
                  gp_datadir=`echo "default(datadir)" | $GP -qf 2>> config.log`
                  AC_MSG_RESULT([$gp_datadir])
                  AC_MSG_CHECKING([comparing GP's and libpari's datadirs])
                  AC_RUN_IFELSE([AC_LANG_PROGRAM([
                      [#include <pari/pari.h>]],
                     [[int t;]
                      [pari_init(5000000, 2);]
                      [t=strcmp(paricfg_datadir,$gp_datadir);]
                      [pari_close()];
                      [return t;]])],
                    [AC_MSG_RESULT([libpari's and GP's datadirs match. Good])
                     AC_MSG_CHECKING([whether pari is configured with pthreads])
                     AC_RUN_IFELSE([AC_LANG_PROGRAM([
                         [#include <pari/pari.h>
                          #include <string.h>]],
                        [[return strcmp(PARI_MT_ENGINE, "pthread") != 0]])],
                       [AC_MSG_RESULT([yes. Good])],
                       [AC_MSG_RESULT([no. Not good])
                        sage_spkg_install_pari=yes],
                       [AC_MSG_RESULT([cross compiling. Assume yes])])
                    ],
                    [AC_MSG_RESULT([libpari's datadir does not match GP's datadir. Not good])
                    sage_spkg_install_pari=yes],
                    [AC_MSG_RESULT([cross compiling. Assume yes])])
                  ], [dnl compared maxver
                   AC_MSG_RESULT([no])
                   sage_spkg_install_pari=yes], [AC_MSG_RESULT([cross compiling. Assume yes])])
                 ], [dnl compared minver
                  AC_MSG_RESULT([no])
                  sage_spkg_install_pari=yes], [AC_MSG_RESULT([cross compiling. Assume yes])])
              AC_LANG_POP()
        ], [sage_spkg_install_pari=yes])
      fi dnl end main PARI test
  ])
  m4_popdef([SAGE_PARI_MAXVER])
  m4_popdef([SAGE_PARI_MINVER])
], [], [], [
    if test x$sage_spkg_install_pari = xyes; then
        AC_SUBST(SAGE_PARI_PREFIX, ['$SAGE_LOCAL'])
    else
       AC_SUBST(SAGE_PARI_PREFIX, [''])
    fi
])
