AC_DEFUN([SAGE_SPKG_ENABLE], [
  m4_pushdef([SPKG_NAME], [$1])
  m4_pushdef([SPKG_TYPE], [$2])
  AC_ARG_ENABLE(SPKG_NAME,
      AS_HELP_STRING([--enable-]SPKG_NAME={no|if_installed|yes},
                     [enable build and use of the $2 package ]SPKG_NAME[ (default: "if_installed")])
AS_HELP_STRING([], [package information: ./sage -info ]SPKG_NAME)
AS_HELP_STRING([--disable-]SPKG_NAME,
                     [disable build and uninstall if previously installed by Sage in PREFIX; same as --enable-]SPKG_NAME[=no]),
      AS_VAR_SET(want_spkg, [$enableval]),
      AS_VAR_SET(want_spkg, [if_installed])
  )

  AS_VAR_SET([is_installed], [no])
  for f in "$SAGE_SPKG_INST/SPKG_NAME"-*; do
      AS_IF([test -r "$f"],
            [AS_VAR_SET([is_installed], [yes])])
  done

  AS_IF([test "$want_spkg" = if_installed],
        [AS_VAR_SET([want_spkg], $is_installed)])

  spkg_line=" \\$(printf '\n    ')SPKG_NAME"
  AS_CASE([$is_installed-$want_spkg],
          [*-yes],  [AS_VAR_APPEND(SAGE_OPTIONAL_INSTALLED_PACKAGES, "$spkg_line")],
          [yes-no], [AS_VAR_APPEND(SAGE_OPTIONAL_CLEANED_PACKAGES, "$spkg_line")])
  AS_VAR_SET([SAGE_ENABLE_]SPKG_NAME, [$want_spkg])
  AC_SUBST([SAGE_ENABLE_]SPKG_NAME)
  m4_popdef([SPKG_TYPE])
  m4_popdef([SPKG_NAME])
])
