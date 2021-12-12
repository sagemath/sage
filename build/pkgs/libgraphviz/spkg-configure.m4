SAGE_SPKG_CONFIGURE([libgraphviz], [
    AC_CHECK_HEADER([graphviz/cgraph.h], [], [sage_spkg_install_libgraphviz=yes])
])
