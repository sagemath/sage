#*****************************************************************************
#       Copyright (C) 2005-2007 William Stein
#       Copyright (C) 2009-2011 David Kirkby
#       Copyright (C) 2012-2016 Jeroen Demeyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

dnl If you are going to update this, please stick the recommended layout
dnl in the autoconf manual - i.e.

dnl First check for programs
dnl Next check for libraries
dnl Next check for header files
dnl Next check for types
dnl Next check for structures
dnl Next check compiler characteristics
dnl Next check for library functions
dnl Next check for system services

dnl Older versions do not support $GFC
AC_PREREQ([2.69])

AC_DEFUN([SAGE_VERSION], m4_esyscmd_s([. src/bin/sage-version.sh && echo $SAGE_VERSION]))
AC_INIT([Sage], SAGE_VERSION, [sage-devel@googlegroups.com])


AC_COPYRIGHT([GPL version 3])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([config])

dnl This is the hierarchy where everything gets installed (SAGE_LOCAL).
dnl By default, it is the subdirectory local, but the user
dnl can override it by supplying --prefix=....
AC_PREFIX_DEFAULT(`pwd -P`/local)

dnl We don't really use automake, but configure needs aclocal and the
dnl automake library files (like config.guess).
AM_INIT_AUTOMAKE([1.9.6 foreign])

# Allow "configure --disable-maintainer-mode" to disable timestamp checking
AM_MAINTAINER_MODE([enable])

dnl Make sure the path to our own m4 macros is always properly set
dnl and doesn't depend on how autoconf is called.
AC_CONFIG_MACRO_DIR([m4])

#---------------------------------------------------------
#
# This is essentially the configure part of the old "install" file.
# Over time, this should be changed to proper autoconf style with
# configure options.
#

########################################################################
# Set various environment variables (needed at configure time)
########################################################################

# Assume current directory is SAGE_ROOT.
SAGE_ROOT=`pwd -P`
AC_SUBST(SAGE_ROOT)

SAGE_LOCAL="$prefix"
# The following line is necessary because the code below
# inspects stuff in the installation prefix at configure time.
# This is nonstandard.
if test "$SAGE_LOCAL" = NONE; then
    SAGE_LOCAL=local
    if test -x "$SAGE_LOCAL"/bin/python3; then
        # Incremental build with an existing installation of python3 spkg
        # in SAGE_LOCAL or venv in SAGE_LOCAL.  Keep old behavior.
        SAGE_VENV_AUTO=no
    else
        SAGE_VENV_AUTO=yes
    fi
else
    SAGE_VENV_AUTO=no
fi
SAGE_SRC="$SAGE_ROOT/src"
SAGE_SPKG_INST="$SAGE_LOCAL/var/lib/sage/installed"

AC_ARG_WITH([sage-venv],
            [AS_HELP_STRING([--with-sage-venv={auto,yes,no,SAGE_VENV}],
                            [put Python packages into an installation hierarchy separate from prefix])],
            [SAGE_VENV="$withval"],
            [SAGE_VENV="auto"])
AC_SUBST([SAGE_VENV])

#---------------------------------------------------------

AC_ARG_ENABLE([build-as-root],
              [AS_HELP_STRING([--enable-build-as-root],
                              [allow building Sage as root (for use in containers)])],
              [enable_build_as_root=$enableval],
              [enable_build_as_root=no])
if test "$enable_build_as_root" != yes; then
   AX_CHECK_ROOT([AC_MSG_ERROR([You cannot build Sage as root, switch to an unprivileged user.  (If building in a container, use --enable-build-as-root.)])], [])
fi

AC_ARG_ENABLE([fat-binary],
              [AS_HELP_STRING([--enable-fat-binary],
                              [build binaries that will run on the widest range of target CPUs (but not relocatable)])],
              [AC_SUBST(SAGE_FAT_BINARY, "yes")],
              [])

AC_ARG_VAR(SAGE_FAT_BINARY, set to "yes" to build binaries that will run on the widest range of target CPUs (but not relocatable))

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug={no|symbols|yes}],
                              [controls debugging support: "no" debugging; debugging "symbols" (default); build debug version ("yes")])],
              [AC_SUBST(SAGE_DEBUG, $enableval)],
              [])

AC_ARG_VAR(SAGE_DEBUG, controls debugging support: "no" debugging; debugging "symbols" (default); build debug version ("yes"))

AC_ARG_ENABLE([editable],
              [AS_HELP_STRING([--enable-editable],
                              [use an editable install of the Sage library])],
              [AC_SUBST([SAGE_EDITABLE], [yes])])

# Check whether we are on a supported platform
AC_CANONICAL_BUILD()
AC_CANONICAL_HOST()

case $host in
*-*-sunos*|*-*-solaris2.[1-9])
AC_MSG_ERROR([[
Sage is not supported on any version of Solaris earlier than 10.
Sage has been tested on the first release of Solaris 10
(03/2005) and works on that. Sage may or may not work with
your version of Solaris.

More information can be found about Sage on Solaris
on the Wiki at http://wiki.sagemath.org/solaris]]);;

*-*-darwin[1-7].*)
AC_MSG_ERROR([[
Sage has never been built on OS X 10.3 (Panther)
or earlier. The oldest version of OS X successfully used
is OS X version 10.4 (Tiger). You might consider updating
your version of OS X if your hardware permits this, but
Apple charges for upgrades of OS X]]);;

*-*-hpux*)
AC_MSG_ERROR([[
You are attempting to build Sage on HP's HP-UX operating system,
which is not a supported platform for Sage yet though
some work has been done on HP-UX. A port does not look to
be particularly difficult. Some information can be
found on the Sage Wiki at http://wiki.sagemath.org/HP-UX

If you would like to help port Sage to HP-UX,
please join the sage-devel discussion list - see
http://groups.google.com/group/sage-devel
The Sage community would also appreciate any patches you submit]]);;

*-*-aix*)
AC_MSG_ERROR([[
You are attempting to build Sage on IBM's AIX operating system,
which is not a supported platform for Sage yet. Things may or
may not work. If you would like to help port Sage to AIX,
please join the sage-devel discussion list - see
http://groups.google.com/group/sage-devel
The Sage community would also appreciate any patches you submit]]);;

*-*-irix*)
AC_MSG_ERROR([[
You are attempting to build Sage on SGI's IRIX operating system,
which is not a supported platform for Sage yet. Things may or
may not work. If you would like to help port Sage to IRIX,
please join the sage-devel discussion list - see
http://groups.google.com/group/sage-devel
The Sage community would also appreciate any patches you submit]]);;

*-*-osf*)
AC_MSG_ERROR([[
You are attempting to build Sage on HP's Tru64 operating system,
which is not a supported platform for Sage yet. Things may or
may not work. If you would like to help port Sage to Tru64,
please join the sage-devel discussion list - see
http://groups.google.com/group/sage-devel
The Sage community would also appreciate any patches you submit]]);;

# The following are all supported platforms.
*-*-freebsd*);;
*-*-linux*);;
*-*-darwin*);;
*-*-solaris*);;
*-*-cygwin*);;

# Wildcard for other unsupported platforms
*)
AC_MSG_ERROR([[
You are attempting to build Sage on $host,
which is not a supported platform for Sage yet. Things may or
may not work. If you would like to help port Sage to $host,
please join the sage-devel discussion list - see
http://groups.google.com/group/sage-devel
The Sage community would also appreciate any patches you submit]]);;
esac


###############################################################################
# Check general programs
###############################################################################

AC_PROG_SED

AC_CHECK_PROG(found_ar, ar, yes, no)
if test x$found_ar != xyes
then
    AC_MSG_NOTICE([Sorry, the 'ar' command must be in the path to build AC_PACKAGE_NAME])
    AC_MSG_NOTICE([On some systems it can be found in /usr/ccs/bin ])
    AC_MSG_NOTICE(['ar' is also part of the GNU 'binutils' package.])
    AC_MSG_ERROR([Exiting, as the archiver 'ar' can not be found.])
fi

AC_CHECK_PROG(found_m4, m4, yes, no)
if test x$found_m4 != xyes
then
    AC_MSG_NOTICE([Sorry, the 'm4' command must be in the path to build AC_PACKAGE_NAME])
    AC_MSG_NOTICE([On some systems it can be found in /usr/ccs/bin])
    AC_MSG_NOTICE([See also http://www.gnu.org/software/m4/])
    AC_MSG_ERROR([Exiting, as the macro processor 'm4' can not be found.])
fi

m4_ifndef([PKG_PROG_PKG_CONFIG], [m4_errprint(
[Error: could not locate the pkg-config autoconf macros. These are
usually located in /usr/share/aclocal/pkg.m4. If your macros are
in a different location, try setting the environment variable
ACLOCAL="aclocal -I/other/macro/dir" before running ./bootstrap
])
dnl Exit autoconf with exit code 16 in this case. This will be
dnl caught by the bootstrap script.
m4_exit(16)])

PKG_PROG_PKG_CONFIG([0.29])

AC_CHECK_PROG(found_ranlib, ranlib, yes, no)
if test x$found_ranlib != xyes
then
    AC_MSG_NOTICE([Sorry, the 'ranlib' command must be in the path to build AC_PACKAGE_NAME])
    AC_MSG_NOTICE([On some systems it can be found in /usr/ccs/bin ])
    AC_MSG_NOTICE(['ranlib' is also part of the GNU 'binutils' package.])
    AC_MSG_ERROR([Exiting, as 'ranlib' can not be found.])
fi

AC_CHECK_PROG(found_strip, strip, yes, no)
if test x$found_strip != xyes
then
    AC_MSG_NOTICE([Sorry, the 'strip' command must be in the path to build AC_PACKAGE_NAME])
    AC_MSG_NOTICE([On some systems 'strip' can be found in /usr/ccs/bin ])
    AC_MSG_NOTICE(['strip' is also part of the GNU 'binutils' package.])
    AC_MSG_ERROR([Exiting, as 'strip' can not be found.])
fi

# Check tar
AC_CACHE_CHECK([for GNU or BSD tar], [ac_cv_path_TAR], [
AC_PATH_PROGS_FEATURE_CHECK(TAR, [tar gtar], [[
ac_version_TAR=`$ac_path_TAR --version 2>&1`
if echo "$ac_version_TAR" | grep >/dev/null GNU; then
    ac_cv_path_TAR=$ac_path_TAR
    if test $ac_prog = tar; then
        ac_path_TAR_found=:
    fi
fi
if echo "$ac_version_TAR" | grep >/dev/null bsdtar; then
    ac_cv_path_TAR=$ac_path_TAR
    if test $ac_prog = tar; then
        ac_path_TAR_found=:
    fi
fi
]],
[AC_MSG_ERROR([could not find either a GNU or BSD version of tar])],
[$PATH:/usr/sfw/bin])
])

command_TAR=`command -v tar 2>/dev/null`
AS_IF([test x$command_TAR != x$ac_cv_path_TAR],
      [AC_MSG_ERROR([[found a good version of tar in $ac_cv_path_TAR, but it's not the first "tar" program in your PATH]])]
)

# Check make (unless MAKE is set)
if test -z "$MAKE"; then
    AC_CACHE_CHECK([for GNU make], [ac_cv_path_MAKE], [
    AC_PATH_PROGS_FEATURE_CHECK(MAKE, [make gmake], [[
    ac_version_MAKE=`$ac_path_MAKE --version 2>&1`
    if echo "$ac_version_MAKE" | grep >/dev/null GNU; then
        ac_cv_path_MAKE=$ac_path_MAKE
        if test $ac_prog = make; then
            ac_path_MAKE_found=:
        fi
    fi
    ]],
    [AC_MSG_ERROR([could not find a GNU version of make])],
    [$PATH:/usr/sfw/bin])
    ])

    command_MAKE=`command -v make 2>/dev/null`
    AS_IF([test x$command_MAKE != x$ac_cv_path_MAKE],
          [AC_MSG_ERROR([[found GNU make in $ac_cv_path_MAKE, but it's not the first "make" program in your PATH]])])
fi

# Check for system python
AC_MSG_CHECKING([for Python])
if SAGE_BOOTSTRAP_PYTHON=$(build/bin/sage-bootstrap-python -c "import sys; print(sys.executable)"); then
    AC_MSG_RESULT([$SAGE_BOOTSTRAP_PYTHON])
else
    AC_MSG_ERROR([You do not have a suitable version of Python installed])
fi

# Check for Latex, the use of which is less important in Sage than
# it used to be, as it was at one time required to build any documentation
# but this is no longer so.
AC_CHECK_PROG(found_latex, latex, yes, no)
if test x$found_latex != xyes
then
    AC_MSG_NOTICE([You do not have 'latex', which is recommended, but not])
    AC_MSG_NOTICE([required. Latex is only really used for building pdf])
    AC_MSG_NOTICE([documents and for %latex mode in the AC_PACKAGE_NAME notebook.])
fi

# Check that perl is available, with version 5.8.0 or later.
# Some packages need perl, however it is not clear whether Sage really
# requires version >= 5.8.0.  The R package *used* to require it, but
# not anymore. -- Jeroen Demeyer
AC_PATH_PROG([PERL],[perl])
AX_PROG_PERL_VERSION([5.8.0],[],[
    AC_MSG_ERROR([Exiting, since AC_PACKAGE_NAME requires perl-5.8.0 or later])
])

###############################################################################
# Check C/C++/Fortran compilers
###############################################################################

# Save compiler flags as configured by the user.
# We have to redo this, because otherwise `AC_PROG_CC` will just overwrite them.
AC_SUBST(CFLAGS, "$CFLAGS")
AC_SUBST(CXXFLAGS, "$CXXFLAGS")
AC_SUBST(FCFLAGS, "$FCFLAGS")
AC_SUBST(F77FLAGS, "$F77FLAGS")
AC_ARG_VAR(CFLAGS, C compiler flags)
AC_ARG_VAR(CXXFLAGS, C compiler flags)
AC_ARG_VAR(FCFLAGS, Fortran compiler flags)
AC_ARG_VAR(F77CFLAGS, Fortran 77 compiler flags)

SAGE_CHECK_CONDA_COMPILERS

AC_PROG_CC()
AC_PROG_CPP()
AC_PROG_CXX()
AC_PROG_FC()

AC_SUBST(CC)
AC_SUBST(FC)

# On darwin, also set the objective C/C++ compilers
# Checking on all platforms doesn't hurt and stops
# configure from sending an error when run on non-darwin.
AC_PROG_OBJC()
AC_PROG_OBJCXX()
AC_SUBST(OBJC)
AC_SUBST(OBJCXX)

AS_IF([test "x$CXX" = x], [AC_MSG_ERROR([a C++ compiler is missing])])


###############################################################################
# Check header files
###############################################################################


# complex.h is one that might not exist on older systems.
AC_LANG(C++)
AC_CHECK_HEADER([complex.h],[],[
    AC_MSG_ERROR([Exiting, since you do not have the 'complex.h' header file.])
])


###############################################################################
# Check libraries
###############################################################################

# First check for something that should be in any maths library (sqrt).
AC_LANG(C++)
AC_CHECK_LIB(m,sqrt,[],[
                       AC_MSG_NOTICE([This system has no maths library installed.])
                       # On AIX this is not installed by default - strange as that might seem.
                       # but is in a fileset bos.adt.libm. However, the fileset bos.adt
                       # includes other things that are probably useful.
                       if test "x`uname`" = 'xAIX'
                       then
                          AC_MSG_NOTICE([On AIX, libm is contained in the bos.adt.libm fileset.])
                          AC_MSG_NOTICE([Actually, we recommend to install the complete bos.adt fileset.])
                          AC_MSG_NOTICE([This needs to be performed by a system administrator.])
                       fi
                       AC_MSG_ERROR([Exiting, since a maths library was not found.])
                       ])

# Check for system services

# Check that we are not building in a directory containing spaces
AS_IF([echo "$ac_pwd" |grep " " >/dev/null],
      AC_MSG_ERROR([the path to the Sage root directory ($ac_pwd) contains a space. Sage will not build correctly in this case])
)

SAGE_CHECK_OSX_SUPPORTED()

###############################################################################
# Collect substitutions for build/make/Makefile.in
###############################################################################

# $(TOOLCHAIN) variable containing prerequisites for the build
SAGE_TOOLCHAIN=gcc
if test "$SAGE_INSTALL_CCACHE" = yes ; then
    SAGE_TOOLCHAIN="$SAGE_TOOLCHAIN ccache"
fi
AC_SUBST([SAGE_TOOLCHAIN])

AC_ARG_ENABLE([experimental-packages],
  [AS_HELP_STRING([--enable-experimental-packages],
    [allow installing experimental packages (default: no = ask for user confirmation for each package)])])
AC_ARG_ENABLE([download-from-upstream-url],
  [AS_HELP_STRING([--enable-download-from-upstream-url],
    [allow downloading packages from their upstream URL if they cannot be found on the Sage mirrors])])

SAGE_SPKG_OPTIONS=""
AS_IF([test "x$enable_experimental_packages" = "xyes"], [
  AS_VAR_APPEND([SAGE_SPKG_OPTIONS], [" -y"])
])
AS_IF([test "x$enable_download_from_upstream_url" = "xyes"], [
  AS_VAR_APPEND([SAGE_SPKG_OPTIONS], [" -o"])
])
AC_SUBST([SAGE_SPKG_OPTIONS])

AC_ARG_ENABLE([cvxopt],
  AS_HELP_STRING([--disable-cvxopt],
                 [disable build of the cvxopt package and its prerequisite suitesparse]), [
    for pkg in cvxopt suitesparse; do
      AS_VAR_SET([SAGE_ENABLE_$pkg], [$enableval])
    done
  ])

AC_ARG_ENABLE([notebook],
  AS_HELP_STRING([--disable-notebook],
                 [disable build of the Jupyter notebook and related packages]), [
    for pkg in notebook nbconvert sagenb_export nbformat terminado send2trash prometheus_client mistune pandocfilters bleach defusedxml jsonschema jupyter_jsmol argon2_cffi; do
      AS_VAR_SET([SAGE_ENABLE_$pkg], [$enableval])
    done
  ])

AS_IF([test "$SAGE_ENABLE_notebook" = no -a "$SAGE_ENABLE_sagelib" = no], [
    for pkg in jupyter_client ipykernel ipython zeromq pyzmq; do
      AS_VAR_SET([SAGE_ENABLE_$pkg], [$enableval])
    done
  ])

AC_ARG_ENABLE([r],
  AS_HELP_STRING([--disable-r],
                 [disable build of the R package and related packages]), [
    for pkg in r rpy2 r_jupyter; do
      AS_VAR_SET([SAGE_ENABLE_$pkg], [$enableval])
    done
  ])

AC_ARG_ENABLE([sagelib],
  AS_HELP_STRING([--disable-sagelib],
                 [disable build of the Sage library and packages depending on it]), [
    for pkg in sagelib sage_docbuild; do
      AS_VAR_SET([SAGE_ENABLE_$pkg], [$enableval])
    done
  ])

dnl Packages that require a special override to use the SPKG
dnl when the system package is not usable.
AS_VAR_SET([sage_use_system_gcc], [force])

SAGE_SPKG_COLLECT()

dnl AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([build/make/Makefile-auto build/make/Makefile])
AC_CONFIG_FILES([src/bin/sage-env-config src/bin/sage-src-env-config build/bin/sage-build-env-config])

AC_CONFIG_FILES([pkgs/sage-conf/sage_conf.py pkgs/sage-conf/setup.cfg])

dnl Create basic directories needed for Sage
AC_CONFIG_COMMANDS(mkdirs,
    [
        for d in "$SAGE_LOGS" "$SAGE_LOCAL" \
            "$SAGE_LOCAL/bin" "$SAGE_LOCAL/etc" \
            "$SAGE_LOCAL/include" "$SAGE_LOCAL/lib" \
            "$SAGE_LOCAL/lib/pkgconfig" \
            "$SAGE_SHARE" "$SAGE_INST"; do
                AC_MSG_NOTICE([creating directory $d])
                mkdir -p "$d" || AC_MSG_ERROR(could not create $d)
        done

        if test -d "$SAGE_LOCAL/lib64"; then
            if test ! -L "$SAGE_LOCAL/lib64"; then
                AC_MSG_NOTICE([[$SAGE_LOCAL/lib64 is not a symlink, see Trac #19782]])
                AC_MSG_ERROR([[Cannot perform incremental update, run "make distclean && make"]])
            fi
        else
            AC_MSG_NOTICE([[creating symbolic link lib64 -> lib]])
            ln -s lib "$SAGE_LOCAL/lib64"
        fi
    ],
    [
        SAGE_LOGS="$SAGE_ROOT/logs/pkgs"
        SAGE_LOCAL="$SAGE_LOCAL"
        SAGE_SHARE="$SAGE_LOCAL/share"
        SAGE_INST="$SAGE_SPKG_INST"
    ])

AC_CONFIG_COMMANDS(links, [
    dnl Create links for the convenience of users
    SYMLINK="${ac_top_build_prefix}prefix"
    AS_IF([test -L "$SYMLINK" -o ! -e "$SYMLINK"], [
        AC_MSG_NOTICE([creating convenience symlink $SYMLINK -> $SAGE_LOCAL])
        ln -sf "$SAGE_LOCAL" "$SYMLINK"
    ], [
        AC_MSG_NOTICE([cannot create convenience symlink $SYMLINK -> $SAGE_LOCAL because the file exists and is not a symlink; this is harmless])
    ])
    SYMLINK="${ac_top_build_prefix}venv"
    AS_IF([test -L "$SYMLINK" -o ! -e "$SYMLINK"], [
        AC_MSG_NOTICE([creating convenience symlink $SYMLINK -> $SAGE_VENV])
        ln -sf "$SAGE_VENV" "$SYMLINK"
    ], [
        AC_MSG_NOTICE([cannot create convenience symlink $SYMLINK -> $SAGE_VENV because the file exists and is not a symlink; this is harmless])
    ])
], [
    SAGE_LOCAL="$SAGE_LOCAL"
    eval SAGE_VENV="$SAGE_VENV"dnl eval so as to evaluate the embedded ${SAGE_LOCAL}
])

AC_OUTPUT()

SAGE_SYSTEM_PACKAGE_NOTICE()

dnl vim:syntax=m4
