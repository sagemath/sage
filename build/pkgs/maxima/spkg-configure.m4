SAGE_SPKG_CONFIGURE([maxima], [
  m4_pushdef([SAGE_MAXIMA_MINVER],["5.45.0"])dnl this version and higher allowed
  SAGE_SPKG_DEPCHECK([ecl], [
    dnl First check for the "maxima" executable in the user's PATH, because
    dnl we still use pexpect to communicate with it in a few places.
    AC_CACHE_CHECK([for Maxima >= $SAGE_MAXIMA_MINVER], [ac_cv_path_MAXIMA], [
        AC_PATH_PROGS_FEATURE_CHECK([MAXIMA], [maxima], [
            maxima_version=`$ac_path_MAXIMA --version 2>&1 | tail -n 1\
                | $SED -n -e 's/Maxima *\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\)/\1/p'`
            AS_IF([test -n "$maxima_version"], [
                AX_COMPARE_VERSION([$maxima_version], [ge], [SAGE_MAXIMA_MINVER], [
                    ac_cv_path_MAXIMA="$ac_path_MAXIMA"
                    ac_path_MAXIMA_found=:
                ])
            ])
        ])
    ])
    SAGE_MAXIMA="$ac_cv_path_MAXIMA"
    AS_IF([test -z "${SAGE_MAXIMA}"], [
      sage_spkg_install_maxima=yes
    ],[
      dnl If we have the executable, check also for the ECL library.
      AC_MSG_CHECKING([if ECL can "require" the maxima module])
      AS_IF([ecl --eval "(require 'maxima)" --eval "(quit)" \
               >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD], [
        AC_MSG_RESULT(yes)
        dnl check also for the Maxima help - needed by Sage
        AC_MSG_CHECKING([if maxima help is working])
	maxima_help_ok=`echo ? ? | ${SAGE_MAXIMA} 2>&1 | grep Details\:`
	AS_IF([test x$maxima_help_ok = x], [AC_MSG_RESULT(yes)], [
               AC_MSG_RESULT(no)
               sage_spkg_install_maxima=yes
        ])
      ], [
        AC_MSG_RESULT(no)
        sage_spkg_install_maxima=yes
      ])
    ])
  ])
  m4_popdef([SAGE_MAXIMA_MINVER])
],[],[],[
  # post-check
  AS_IF([test x$sage_spkg_install_maxima = xyes], [
    dnl Leaving this variable empty will tell sagelib to load
    dnl the maxima library (within ECL) by name instead of by
    dnl absolute path.
    SAGE_MAXIMA='${prefix}'/bin/maxima
    SAGE_MAXIMA_FAS='${prefix}'/lib/ecl/maxima.fas
  ])

  AC_SUBST(SAGE_MAXIMA, "${SAGE_MAXIMA}")
  AC_SUBST(SAGE_MAXIMA_FAS, "${SAGE_MAXIMA_FAS}")
])
