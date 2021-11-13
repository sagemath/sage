SAGE_SPKG_CONFIGURE([maxima], [
  SAGE_SPKG_DEPCHECK([ecl], [
    dnl First check for the "maxima" executable in the user's PATH, because
    dnl we still use pexpect to communicate with it in a few places.
    AC_PATH_PROG([SAGE_MAXIMA], [maxima])
    AS_IF([test -z "${SAGE_MAXIMA}"], [
      sage_spkg_install_maxima=yes
    ],[
      dnl If we have the executable, check also for the ECL library.
      AC_MSG_CHECKING([if ECL can "require" the maxima module])
      AS_IF([ecl --eval "(require 'maxima)" --eval "(quit)" \
               >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD], [
        AC_MSG_RESULT(yes)
      ], [
        AC_MSG_RESULT(no)
        sage_spkg_install_maxima=yes
      ])
    ])
  ])
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
