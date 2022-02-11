SAGE_SPKG_CONFIGURE([graphviz], [
    dnl We check all executables that are tested by sage.features.graphviz
    AC_CHECK_PROGS([DOT], [dot])
    AS_IF([test x$DOT = x], [sage_spkg_install_graphviz=yes])
    AC_CHECK_PROGS([NEATO], [neato])
    AS_IF([test x$NEATO = x], [sage_spkg_install_graphviz=yes])
    AC_CHECK_PROGS([TWOPI], [twopi])
    AS_IF([test x$TWOPI = x], [sage_spkg_install_graphviz=yes])
])
