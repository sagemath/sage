AC_INIT([ATLAS],[3.10.1])
AM_INIT_AUTOMAKE(foreign)

AC_CANONICAL_HOST

LT_INIT

AC_PROG_CC
AC_PROG_F77

AC_CONFIG_FILES([Makefile])
AC_CONFIG_MACRO_DIR([m4])

AC_CHECK_LIB(pthread, pthread_create, [PTHREAD_LIB="-lpthread"])
AC_SUBST(PTHREAD_LIB)

AC_ARG_ENABLE([static],
    [AS_HELP_STRING([--enable-static],
                    [also install static libtool library (default: no)])])
AS_IF([test "x$enable_static" = "xyes"],
    [libtool_type=""],
    [libtool_type=-shared])
AC_SUBST(LIBTOOL_TYPE, [$libtool_type])

case $host in
     *-*-cygwin*)
        LTPREFIX=cyg
        LTSUFFIX=.dll
        LTCLINK="\$(CC) -L. -shared"
        LTF77LINK="\$(F77) -L. -shared"
        LTLINKFLAGS="-Wl,--out-implib=lib\$@.dll.a -Wl,--export-all-symbols -Wl,--enable-auto-import"
        LTINSTALL="\$(CP) \$(LTPREFIX)\$\${module}\$(LTSUFFIX) \$(bindir)/\$(LTPREFIX)\$\${module}\$(LTSUFFIX); \$(CP) lib\$\${module}.dll.a \$(libdir)/lib\$\${module}.dll.a"
        ;;
     *)
        LTPREFIX=lib
        LTSUFFIX=.la
        LTCLINK="\$(LIBTOOL) --tag=CC --mode=link \$(CC) \$(LIBTOOL_FLAG)"
        LTF77LINK="\$(LIBTOOL) --tag=F77 --mode=link \$(F77) \$(LIBTOOL_FLAG)"
        LTLINKFLAGS="-rpath \$(libdir) -version-info \$(SO_VERSION)"
        LTINSTALL="\$(LIBTOOL) --mode=install \$(INSTALL) -c \$(LTPREFIX)\$\${module}\$(LTSUFFIX) \$(libdir)/\$(LTPREFIX)\$\${module}\$(LTSUFFIX)"
        ;;
esac

AC_SUBST(LTPREFIX)
AC_SUBST(LTSUFFIX)
AC_SUBST(LTCLINK)
AC_SUBST(LTF77LINK)
AC_SUBST(LTLINKFLAGS)
AC_SUBST(LTINSTALL)

AC_OUTPUT
