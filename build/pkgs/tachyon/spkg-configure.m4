SAGE_SPKG_CONFIGURE([tachyon], [
        AC_PATH_PROG([TACHYON], [tachyon])
        AS_IF([test x$TACHYON = x], [
           AC_MSG_NOTICE([tachyon not found. Installing tachyon])
           sage_spkg_install_tachyon=yes])
])
