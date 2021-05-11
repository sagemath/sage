SAGE_SPKG_CONFIGURE([coxeter3], [
    AC_LANG_PUSH(C++)
    AC_MSG_CHECKING([for library coxeter3])
    SAVE_LIBS="$LIBS"
    LIBS="$LIBS -lcoxeter3"
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[
            #include <coxeter/sage.h>
            #include <coxeter/interactive.h>
          ]], [[
           coxeter::CoxGroup *g = interactive::coxeterGroup("B", 2);
          ]])
      ], [
        AC_MSG_RESULT([yes])
      ], [
        AC_MSG_RESULT([no])
        sage_spkg_install_coxeter3=yes
      ])
    LIBS="$SAVE_LIBS"
    AC_LANG_POP(C++)
])
