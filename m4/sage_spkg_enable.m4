AC_DEFUN([SAGE_SPKG_ENABLE], [dnl
  m4_pushdef([SPKG_NAME], [$1])dnl
  m4_pushdef([SPKG_TYPE], [$2])dnl
  m4_if(SPKG_TYPE, [standard], [dnl
    dnl standard packages; help message is deliberately very brief,
    dnl as this is for advanced users only
    AC_ARG_ENABLE(SPKG_NAME,
      AS_HELP_STRING([--disable-]SPKG_NAME,
                     [disable $2 package ]SPKG_NAME),
      AS_VAR_SET([SAGE_ENABLE_]SPKG_NAME, [$enableval])
    )
  ], [dnl
    dnl optional/experimental packages
    AC_ARG_ENABLE(SPKG_NAME,
      AS_HELP_STRING([--enable-]SPKG_NAME={no|if_installed (default)|yes},
                     [enable build and use of the SPKG_TYPE package $3], [26], [100])
AS_HELP_STRING([], [* package info: ./sage -info SPKG_NAME])
AS_HELP_STRING([--disable-]SPKG_NAME,
                     [disable build and uninstall if previously installed by Sage in PREFIX; same as --enable-]SPKG_NAME[=no], [26], [100]),
      AS_VAR_SET([SAGE_ENABLE_]SPKG_NAME, [$enableval])
    )
  ])dnl
  m4_popdef([SPKG_TYPE])dnl
  m4_popdef([SPKG_NAME])dnl
])
