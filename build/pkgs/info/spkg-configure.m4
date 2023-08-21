SAGE_SPKG_CONFIGURE([info], [
  AC_PATH_PROG(INFO, info)
  AS_IF([test -z "${INFO}"], [sage_spkg_install_info=yes
    ], [
     dnl very old makeinfo are not texi2any, newer are symlinks to texi2any
     AC_PATH_PROG(TEXI2ANY, texi2any)
     AS_IF([test -z "${TEXI2ANY}"], [sage_spkg_install_info=yes
       ], [
        AS_IF([makeinfo -c foo 2>&1 | grep -q invalid], [
         dnl makeinfo found, but too old, and  does not support all options that ecl likes to use
         sage_spkg_install_info=yes])
        rm -f stdin.info
       ])
    ])
])
