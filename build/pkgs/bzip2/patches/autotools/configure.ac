#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.57])
AC_INIT([bzip2], [1.0.6], [Julian Seward <jseward@bzip.org>])
BZIP2_LT_CURRENT=1
BZIP2_LT_REVISION=6
BZIP2_LT_AGE=0
AC_CONFIG_SRCDIR([bzlib.h])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
LT_INIT([win32-dll])
PKG_PROG_PKG_CONFIG

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Check for system features.
AC_SYS_LARGEFILE

AC_MSG_CHECKING([whether compiler understands -Wall])
save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -Wall"
AC_TRY_COMPILE([], [], [
	AC_MSG_RESULT([yes])
], [
	AC_MSG_RESULT([no])
	CFLAGS="$save_CFLAGS"
])

AC_MSG_CHECKING([whether compiler understands -Winline])
save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -Winline"
AC_TRY_COMPILE([], [], [
	AC_MSG_RESULT([yes])
], [
	AC_MSG_RESULT([no])
	CFLAGS="$save_CFLAGS"
])

AC_MSG_CHECKING([whether compiler understands -fno-strength-reduce])
save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -fno-strength-reduce"
AC_TRY_COMPILE([], [], [
	AC_MSG_RESULT([yes])
], [
	AC_MSG_RESULT([no])
	CFLAGS="$save_CFLAGS"
])

# Checks for library functions.

# Write the output.
AC_SUBST([BZIP2_LT_CURRENT])
AC_SUBST([BZIP2_LT_REVISION])
AC_SUBST([BZIP2_LT_AGE])
AC_CONFIG_FILES([Makefile bzip2.pc])
AC_OUTPUT
