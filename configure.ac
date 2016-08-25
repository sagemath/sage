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

dnl We don't really use automake, but configure needs aclocal and the
dnl automake library files (like config.guess).
AM_INIT_AUTOMAKE([1.9.6 foreign])

# Allow "configure --disable-maintainer-mode" to disable timestamp checking
AM_MAINTAINER_MODE([enable])

dnl Make sure the path to our own m4 macros is always properly set 
dnl and doesn't depend on how autoconf is called.
AC_CONFIG_MACRO_DIR([m4])

#---------------------------------------------------------
# We need to run this configure script with bash
if test -z "$BASH_VERSION$CONFIG_SHELL"
then
    CONFIG_SHELL=bash
    export CONFIG_SHELL
    if $CONFIG_SHELL -c "exit 0"
    then
        exec $CONFIG_SHELL $0 "$@"
    else
        AC_MSG_NOTICE([The 'bash' shell is needed to build AC_PACKAGE_NAME])
        AC_MSG_NOTICE([All modern systems will have the 'bash' shell installed somewhere])
        if test -d /opt/OpenSource/bin
        then
           AC_MSG_NOTICE([On HP-UX you may try adding /opt/OpenSource/bin to your path])
        fi
        if test -d /opt/pware/bin
        then
           AC_MSG_NOTICE([On AIX you may try adding /opt/pware/bin to your path])
        fi
        AC_MSG_ERROR(['bash' not found])
    fi
fi

#---------------------------------------------------------
#
# This is essentially the configure part of the old "install" file.
# Over time, this should be changed to proper autoconf style with
# configure options.
#

########################################################################
# Set various environment variables
########################################################################

# Assume current directory is SAGE_ROOT
SAGE_ROOT=`pwd -P`
SAGE_SRC="$SAGE_ROOT/src"
SAGE_LOCAL="$SAGE_ROOT/local"
SAGE_SHARE="$SAGE_LOCAL/share"
SAGE_EXTCODE="$SAGE_SHARE/sage/ext"
SAGE_SPKG_INST="$SAGE_LOCAL/var/lib/sage/installed"

#---------------------------------------------------------
AX_CHECK_ROOT([AC_MSG_ERROR([You cannot build Sage as root, switch to an unpriviledged user])], [])

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

*-*-freebsd*)
AC_MSG_ERROR([[
You are attempting to build Sage on the FreeBSD operating system,
which is not a supported platform for Sage yet, though
developers are working on adding FreeBSD support. Things may or
may not work. If you would like to help port Sage to FreeBSD,
please join the sage-devel discussion list - see
http://groups.google.com/group/sage-devel
The Sage community would also appreciate any patches you submit]]);;

# The following are all supported platforms.
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

# Check for Latex, the use of which is less important in Sage than
# it used to be, as it was at one time required to build any documentation
# but this is no longer so.
AC_CHECK_PROG(found_latex, latex, yes, no)
if test x$found_latex != xyes
then
    AC_MSG_WARN([You do not have 'latex', which is recommended, but not])
    AC_MSG_WARN([required. Latex is only really used for building pdf])
    AC_MSG_WARN([documents and for %latex mode in the AC_PACKAGE_NAME notebook.])
fi

# Check that perl is available, with version 5.8.0 or later.
# Some packages need perl, however it is not clear whether Sage really
# requires version >= 5.8.0.  The R package *used* to require it, but
# not anymore. -- Jeroen Demeyer
AC_PATH_PROG([PERL],[perl])
AX_PROG_PERL_VERSION([5.8.0],[],[
    AC_MSG_ERROR([Exiting, since AC_PACKAGE_NAME requires perl-5.8.0 or later])
])

# Check whether git works by executing "git --version"
need_to_install_git=no
AC_CACHE_CHECK([for git], [ac_cv_path_GIT], [
    AC_PATH_PROGS_FEATURE_CHECK([GIT], [git],
    [${ac_path_GIT} --version >/dev/null 2>/dev/null && ac_cv_path_GIT=${ac_path_GIT}],
    [need_to_install_git=yes; ac_cv_path_GIT=no])])


###############################################################################
# Check C/C++/Fortran compilers
###############################################################################

dnl Usage: SAGE_SHOULD_INSTALL_GCC(reason)
dnl
dnl Use this macro to indicate that we SHOULD install GCC.
dnl In this case, GCC will be installed unless SAGE_INSTALL_GCC=no.
dnl In the latter case, a warning is given.
AC_DEFUN([SAGE_SHOULD_INSTALL_GCC], [
    if test x$SAGE_INSTALL_GCC = xexists; then
        true  # Do nothing if already installed
    elif test x$SAGE_INSTALL_GCC = xno; then
        AC_MSG_WARN([$1])
    else
        AC_MSG_NOTICE([Installing GCC because $1])
        need_to_install_gcc=yes
    fi
])

dnl Usage: SAGE_MUST_INSTALL_GCC(reason)
dnl
dnl Use this macro to indicate that we MUST install GCC.
dnl In this case, it is an error if SAGE_INSTALL_GCC=no.
AC_DEFUN([SAGE_MUST_INSTALL_GCC], [
    if test x$SAGE_INSTALL_GCC = xexists; then
        true  # Do nothing if already installed
    elif test x$SAGE_INSTALL_GCC = xno; then
        AC_MSG_ERROR([SAGE_INSTALL_GCC is set to 'no', but $1])
    else
        AC_MSG_NOTICE([Installing GCC because $1])
        need_to_install_gcc=yes
    fi
])


# By default, do not install GCC
need_to_install_gcc=no

if test -f "$SAGE_LOCAL/bin/gcc"; then
    # Special value for SAGE_INSTALL_GCC if GCC is already installed
    SAGE_INSTALL_GCC=exists
elif test -n "$SAGE_INSTALL_GCC"; then
    # Check the value of the environment variable SAGE_INSTALL_GCC
    case "$SAGE_INSTALL_GCC" in
        yes)
            SAGE_MUST_INSTALL_GCC([SAGE_INSTALL_GCC is set to 'yes']);;
        no)
            true;;
        *)
            AC_MSG_ERROR([SAGE_INSTALL_GCC should be set to 'yes' or 'no'. You can also leave it unset to install GCC when needed]);;
    esac
fi


AC_LANG(C)
AC_PROG_CC()
AC_PROG_CPP()

AC_LANG(C++)
AC_PROG_CXX()

AC_LANG(Fortran)
AC_PROG_FC()

if test "x$CXX" = x
then
    AC_MSG_ERROR([a C++ compiler is missing])
fi


###############################################################################
# Check header files
###############################################################################

# complex.h is one that might not exist on older systems.
AC_LANG(C++)
AC_CHECK_HEADER([complex.h],[],[
    AC_MSG_ERROR([Exiting, since you do not have the 'complex.h' header file.])
])


###############################################################################
# Check types/structures
###############################################################################

# None needed


###############################################################################
# Check compiler characteristics
###############################################################################

AX_CXX_COMPILE_STDCXX_11([], optional)
if test $HAVE_CXX11 != 1; then
    SAGE_MUST_INSTALL_GCC([your C++ compiler does not support C++11])
fi

AC_LANG(C)
if test -z "$CC"; then
    SAGE_MUST_INSTALL_GCC([a C compiler is missing])
else
    # Check that we can compile C99 code
    AC_PROG_CC_C99()
    if test "x$ac_cv_prog_cc_c99" = xno
    then
        SAGE_MUST_INSTALL_GCC([your C compiler cannot compile C99 code])
    fi
fi

# Check that the Fortran compiler accepts free-format source code
# (as opposed to the older fixed-format style from Fortran 77).
# This helps verify the compiler works too, so if some idiot
# sets FC to /usr/bin/ls, we will at least know it's
# not a working Fortran compiler.
AC_LANG(Fortran)
if test -z "$FC"; then
    SAGE_MUST_INSTALL_GCC([a Fortran compiler is missing])
else
   # see http://www.gnu.org/software/hello/manual/autoconf/Fortran-Compiler.html
   AC_FC_FREEFORM([],
   [
   AC_MSG_NOTICE([Your Fortran compiler does not accept free-format source code])
   AC_MSG_NOTICE([which means the compiler is either seriously broken, or])
   AC_MSG_NOTICE([is too old to build Sage.])
   SAGE_MUST_INSTALL_GCC([the Fortran compiler is not suitable])
   ])
fi


# Check compiler versions

if test x$GXX != xyes; then
    SAGE_SHOULD_INSTALL_GCC([your C++ compiler isn't GCC (GNU C++)])
elif test $need_to_install_gcc = yes; then
    # If we're installing GCC anyway, skip the rest of these version
    # checks.
    true
elif test x$GCC != xyes; then
    SAGE_SHOULD_INSTALL_GCC([your C compiler isn't GCC (GNU C)])
elif test x$GFC != xyes; then
    SAGE_SHOULD_INSTALL_GCC([your Fortran compiler isn't GCC (GNU Fortran)])
else
    # Since need_to_install_gcc is "no", we know that
    # at least C, C++ and Fortran compilers are available.
    # We also know that all compilers are GCC.

    # Find out the compiler versions:
    AX_GCC_VERSION()
    AX_GXX_VERSION()

    # Add the .0 because Debian/Ubuntu gives version numbers like
    # 4.6 instead of 4.6.4 (Trac #18885)
    case "$GXX_VERSION.0" in
        [[0-3]].*|4.[[0-7]].*)
            # Install our own GCC if the system-provided one is older than gcc-4.8.
            SAGE_SHOULD_INSTALL_GCC([you have $CXX version $GXX_VERSION, which is quite old]);;
    esac

    # The following tests check that the version of the compilers
    # are all the same.
    if test "$GCC_VERSION" != "$GXX_VERSION"; then
        SAGE_SHOULD_INSTALL_GCC([$CC ($GCC_VERSION) and $CXX ($GXX_VERSION) are not the same version])
    fi

    # In the paragraph below, 'gfortran' is used to indicate the GNU Fortran
    # compiler, though it might be called something else.

    # It's not easily possible to determine the Fortran version, as
    # gfortran -dumpversion did not until GCC 4.5 return just the
    # the version number, but the same as gfortran --version
    # for example:

    # drkirkby@hawk:~$ gcc -dumpversion
    # 4.3.4

    # drkirkby@hawk:~$ g++ -dumpversion
    # 4.3.4

    # drkirkby@hawk:~$ gfortran -dumpversion
    # GNU Fortran (GCC) 4.3.4
    # Copyright (C) 2008 Free Software Foundation, Inc.
    # GNU Fortran comes with NO WARRANTY, to the extent permitted by law.
    # You may redistribute copies of GNU Fortran
    # under the terms of the GNU General Public License.
    # For more information about these matters, see the file named COPYING

    # This behaviour is fixed in the gcc 4.5 branch. Since we need to
    # support older versions of the compiler, we can't change this.

    # But I would expect that the version will be on the output
    # of the compiler followed by -dumpversion (e.g. fortran -dumpversion)

    # So we grep for the known gcc version on the output of gfortran -dumpversion.
    # and hope we find the same string. If so, they are almost certainly
    # the same version.
    fortran_version_string="`$FC -dumpversion | grep $GCC_VERSION 2>&1`"

    if test "x$fortran_version_string" = x
    then
        SAGE_SHOULD_INSTALL_GCC([the Fortran compiler is not the same version as the C compiler])
    fi
fi


# Check that the assembler and linker used by $CXX match $AS and $LD.
# See http://trac.sagemath.org/sage_trac/ticket/14296
if test -n "$AS"; then
    CXX_as=`$CXX -print-file-name=as 2>/dev/null`
    CXX_as=`command -v $CXX_as 2>/dev/null`
    cmd_AS=`command -v $AS`

    if test "$CXX_as" != "" -a "$CXX_as" != "$cmd_AS"; then
        SAGE_SHOULD_INSTALL_GCC([there is a mismatch of assemblers])
        AC_MSG_NOTICE([  $CXX uses $CXX_as])
        AC_MSG_NOTICE([  \$AS equal to $AS])
    fi
fi
if test -n "$LD"; then
    CXX_ld=`$CXX -print-file-name=ld 2>/dev/null`
    CXX_ld=`command -v $CXX_ld 2>/dev/null`
    cmd_LD=`command -v $LD`

    if test "$CXX_ld" != "" -a "$CXX_ld" != "$cmd_LD"; then
        SAGE_SHOULD_INSTALL_GCC([there is a mismatch of linkers])
        AC_MSG_NOTICE([  $CXX uses $CXX_ld])
        AC_MSG_NOTICE([  \$LD equal to $LD])
    fi
fi


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


if test x`uname` = xDarwin; then
[
    # Warning: xcodebuild does not seem to be maintained in Xcode 4.3
    # or later, so do not rely on the variable XCODE_VERS with OS X
    # 10.7 or later.
    XCODE_VERS=`xcodebuild -version 2> /dev/null | grep Xcode | sed -e 's/[A-Za-z ]//g'`
    if [ -z $XCODE_VERS ]; then
        XCODE_VERS="2"
    fi
    XCODE_VERS_MAJOR=`echo $XCODE_VERS | cut '-d.' -f1`
    DARWIN_VERSION=`uname -r | cut '-d.' -f1`
    echo "***************************************************"
    echo "***************************************************"
    if [ $DARWIN_VERSION -gt 10 ]; then
        echo "You are using OS X Lion (or later)."
        echo "You are strongly advised to install Apple's latest Xcode"
        echo "unless you already have it. You can install this using"
        echo "the App Store. Also, make sure you install Xcode's"
        echo "Command Line Tools -- see Sage's README.txt."
    elif [ $XCODE_VERS_MAJOR -gt 2 ]; then
        echo "You are using Xcode version $XCODE_VERS."
        echo "You are strongly advised to install Apple's latest Xcode"
        echo "unless you already have it. You can download this from"
        echo "http://developer.apple.com/downloads/."
        echo "If using Xcode 4.3 or later, make sure you install Xcode's"
        echo "Command Line Tools -- see Sage's README.txt."
    else
        echo "You are using Xcode version 1 or 2"
        echo "WARNING: You are strongly advised to install the"
        echo "latest version of Apple's Xcode for your platform,"
        echo "unless you already have it."
        if [ $DARWIN_VERSION -eq 10 ]; then
           echo "Probably you need Xcode 3.2.6"
        elif [ $DARWIN_VERSION -eq 9 ]; then
           echo "Probably you need Xcode 3.1.4"
        elif [ $DARWIN_VERSION -lt 9 ]; then
           echo "Probably you need Xcode 2.5"
        fi
    fi
]

###########################################################################
# (OS X only)
# Sage will probably not build at all if either Fink or MacPorts can be
# found, and the error messages can be extremely confusing.  Even if it does
# build, the product will probably be wrong.  This runs a basic check to
# find them. Once the Sage build process is perfected, this won't be necessary.
# dphilp 15/9/2008
###########################################################################
    PORTS_PATH=`which port`
    if test -f "$PORTS_PATH"; then
AC_MSG_ERROR([["found MacPorts in $PORTS_PATH. Either:
(1) rename /opt/local and /sw, or
(2) change PATH and DYLD_LIBRARY_PATH
(Once Sage is built, you can restore them.)]])
    fi

    FINK_PATH=`which fink`
    if test -f "$FINK_PATH"; then
AC_MSG_ERROR([["found Fink in $FINK_PATH. Either:
(1) rename /opt/local and /sw, or
(2) change PATH and DYLD_LIBRARY_PATH
(Once Sage is built, you can restore them.)]])
    fi
fi


###############################################################################
# Create $SAGE_ROOT/build/make/Makefile starting from build/make/deps
###############################################################################

# Use file descriptor 7 since make uses 3 and 4 and configure uses 5 and 6
exec 7>build/make/Makefile

cat >&7 <<EOF
#==============================================================================
# This file has been automatically generated by
#   $SAGE_ROOT/configure
# You should not edit it by hand
#==============================================================================

# Always use bash for make rules
SHELL = `command -v bash`

# Directory to keep track of which packages are installed
INST = `echo "$SAGE_SPKG_INST" | sed 's/ /\\\\ /g'`

# Aliases for optional packages selected at configure time
EOF


# Python version
if test "$SAGE_PYTHON3" = yes; then
    echo >&7 'PYTHON = $(inst_python3)'
else
    echo >&7 'PYTHON = $(inst_python2)'
fi

# Sage MP library
AC_MSG_CHECKING([multiprecision library])
AC_ARG_WITH([mp],
[AS_HELP_STRING([--with-mp=mpir],
    [use MPIR as multiprecision library (default)])]
[AS_HELP_STRING([--with-mp=gmp],
    [use GMP as multiprecision library])])

case "$with_mp" in
    MPIR|mpir) with_mp=MPIR;;
    GMP|gmp) with_mp=GMP;;
    "") with_mp=${SAGE_MP_LIBRARY:-MPIR};;
    *)
        AC_MSG_ERROR([allowed values for --with-mp are MPIR and GMP]) ;;
esac

AC_MSG_RESULT([$with_mp])

if test $with_mp = MPIR; then
    echo >&7 'MP_LIBRARY = $(inst_mpir)'
else
    echo >&7 'MP_LIBRARY = $(inst_gmp)'
fi

# BLAS library
AC_MSG_CHECKING([BLAS library])
AC_ARG_WITH([blas],
[AS_HELP_STRING([--with-blas=atlas],
    [use ATLAS as BLAS library (default)])]
[AS_HELP_STRING([--with-blas=openblas],
    [use OpenBLAS as BLAS library])],,
    [with_blas=atlas]  # default
)

case "$with_blas" in
    atlas) true;;
    openblas) need_to_install_atlas=no;;
    *) AC_MSG_ERROR([allowed values for --with-blas are 'atlas' and 'openblas']);;
esac

AC_MSG_RESULT([$with_blas])

echo >&7 "BLAS = \$(inst_${with_blas})"

# $(TOOLCHAIN) variable containing prerequisites for the build
AS_ECHO_N(>&7 ['TOOLCHAIN = $(inst_gcc)'])
if test "$SAGE_INSTALL_CCACHE" = yes ; then
    AS_ECHO_N(>&7 [' $(inst_ccache)'])
fi
echo >&7

AC_MSG_CHECKING([package versions])
AC_MSG_RESULT([])

# Usage: newest_version $pkg
# Print version number of latest package $pkg
newest_version() {
    PKG=$1
    if test -f "$SAGE_ROOT/build/pkgs/$PKG/package-version.txt" ; then
        AS_ECHO_N(["$PKG-"])
        cat "$SAGE_ROOT/build/pkgs/$PKG/package-version.txt"
    else
        echo "$PKG"
    fi
}

# Outputs the list of packages, filtered by 'type', e.g.:
#
#     filtered_packages_list base
#     filtered_packages_list standard
#     filtered_packages_list optional
#     filtered_packages_list experimental
#
# Or, if you want all packages:
#
#     filtered_packages_list all
#
# Or, if you want all packages which should appear in the source tarball:
#
#     filtered_packages_list sdist
#
# The output consists of triples
# PKG_NAME PKG_VERSION PKG_VAR
#

changequote(<,>)
filtered_packages_list() {
    # for each package in pkgs/
    for DIR in $SAGE_ROOT/build/pkgs/*; do
        test -d "$DIR" || continue

        PKG_TYPE_FILE="$DIR/type"
        if [ -f "$PKG_TYPE_FILE" ]; then
            PKG_TYPE=`cat $PKG_TYPE_FILE`
        else
            # exit won't necessarily exit 'configure', so signal an error this way, too:
            echo INVALID
            changequote([,])
            AC_MSG_ERROR(["$PKG_TYPE_FILE" is missing.])
            changequote(<,>)
        fi

        # Check consistency of 'DIR/type' file
        if [ "$PKG_TYPE" != "base" ] && \
           [ "$PKG_TYPE" != "standard" ] && \
           [ "$PKG_TYPE" != "optional" ] && \
           [ "$PKG_TYPE" != "experimental" ] && \
           [ "$PKG_TYPE" != "script" ] && \
           [ "$PKG_TYPE" != "pip" ]; then
            # exit won't necessarily exit 'configure', so signal an error this way, too:
            echo INVALID
            changequote([,])
            AC_MSG_ERROR([The content of "$PKG_TYPE_FILE" must be 'base', 'standard', 'optional', 'experimental', 'script', or 'pip'])
            changequote(<,>)
        fi

        PKG_NAME=$(basename $DIR)

        # Filter
        filter=no
        if [ "$1" = all ]; then
            filter=yes
        elif [ "$1" = "$PKG_TYPE" ]; then
            filter=yes
        elif [ "$1" = sdist ]; then
            # sdist packages are all standard packages, together with
            # mpir and python2.
            if [ "$PKG_TYPE" = standard ]; then
                filter=yes
            elif [ "$PKG_NAME" = mpir ]; then
                filter=yes
            elif [ "$PKG_NAME" = python2 ]; then
                filter=yes
            fi
        fi
        if [ $filter = yes ]; then
            PKG_VAR="$(echo $PKG_NAME | sed 's/^/inst_/')"
            PKG_VERSION=$(newest_version $PKG_NAME)
            echo "$PKG_NAME $PKG_VERSION $PKG_VAR"
        fi
    done
}
changequote([,])

echo >&7
echo >&7 '# Files to track installation of packages'
filtered_packages_list all | while read PKG_NAME PKG_VERSION PKG_VAR; do
    if test "$PKG_NAME" = INVALID; then
        # filtered_packages_list signals an error (already reported there)
        exit 1 # this doesn't leave 'configure' (yet)
    elif test "$PKG_NAME" != "$PKG_VERSION"; then
        # If $need_to_install_{PKG_NAME} is set to no, then set $PKG_VAR
        # to some dummy file to skip the installation. Note that an
        # explicit "./sage -i PKG_NAME" will still install the package.
        need_to_install="need_to_install_${PKG_NAME}"
        if test "${!need_to_install}" != no ; then
            echo >&7 "$PKG_VAR = \$(INST)/$PKG_VERSION"
            AC_MSG_RESULT([    $PKG_VERSION])
        else
            echo >&7 "$PKG_VAR = \$(INST)/.dummy"
            AC_MSG_RESULT([    $PKG_VERSION not installed (configure check)])
        fi
    fi
done || exit 1

echo >&7
echo >&7

[
echo >&7 '# All standard packages'
echo >&7 'STANDARD_PACKAGES = \'
filtered_packages_list standard | while read PKG_NAME PKG_VERSION PKG_VAR; do
    if test "$PKG_NAME" = INVALID; then
        # filtered_packages_list signals an error (already reported there)
        exit 1 # this doesn't leave 'configure' (yet)
    fi
    echo >&7 "    \$($PKG_VAR) \\"
done || exit 1
echo >&7
echo >&7

echo >&7 '# All optional installed packages (triggers the auto-update)'
echo >&7 'OPTIONAL_INSTALLED_PACKAGES = \'
filtered_packages_list optional | while read PKG_NAME PKG_VERSION PKG_VAR; do
    if test "$PKG_NAME" = INVALID; then
        # filtered_packages_list signals an error (already reported there)
        exit 1 # this doesn't leave 'configure' (yet)
    fi
    if [ -f $SAGE_SPKG_INST/$PKG_NAME-* ]; then
        echo >&7 "    \$($PKG_VAR) \\"
    fi;
done || exit 1
echo >&7
echo >&7

echo >&7 '# All packages which should be downloaded'
echo >&7 'SDIST_PACKAGES = \'
filtered_packages_list sdist | while read PKG_NAME PKG_VERSION PKG_VAR; do
    if test "$PKG_NAME" = INVALID; then
        # filtered_packages_list signals an error (already reported there)
        exit 1 # this doesn't leave 'configure' (yet)
    fi
    echo >&7 "    \$(INST)/$PKG_VERSION \\"
done || exit 1
echo >&7
echo >&7


echo >&7 'SCRIPT_SOURCES = \'
for file in "$SAGE_SRC/bin/"*; do
    echo >&7 "                 \$(SAGE_SRC)${file#$SAGE_SRC} \\"
done
echo >&7
echo >&7 'SCRIPTS = \'
for file in "$SAGE_SRC/bin/"*; do
    echo >&7 "                 \$(SAGE_LOCAL)${file#$SAGE_SRC} \\"
done
echo >&7
echo >&7 'EXTCODE_SOURCES = \'
for file in `find "$SAGE_SRC"/ext -type f`; do
    echo >&7 "                  \$(SAGE_SRC)${file#$SAGE_SRC} \\"
done
echo >&7
echo >&7 'EXTCODE = \'
for file in `find "$SAGE_SRC"/ext -type f`; do
    echo >&7 "                  \$(SAGE_EXTCODE)${file#$SAGE_SRC/ext} \\"
done
echo >&7

# Copy build/make/deps
cat >&7 <<EOF

# Dummy target for packages which are not installed
\$(INST)/.dummy:
	touch \$@

#==============================================================================
# What follows now is a copy of
#   $SAGE_ROOT/build/make/deps
#==============================================================================

EOF

cat "$SAGE_ROOT/build/make/deps" >&7
# Copy build/make/deps

cat >&7 <<EOF

#==============================================================================
# Rules generated from pkgs/<package>/dependencies files
#==============================================================================

EOF

# Add a Makefile target corresponding to a given package
filtered_packages_list all | {
while read PKG_NAME PKG_VERSION PKG_VAR; do
    if test "$PKG_NAME" = INVALID; then
        # filtered_packages_list signals an error (already reported there)
        exit 1 # this doesn't leave 'configure' (yet)
    fi
    DEP_FILE="$SAGE_ROOT/build/pkgs/$PKG_NAME/dependencies"
    TYPE=`cat "$SAGE_ROOT/build/pkgs/$PKG_NAME/type"`
    if [ -f "$DEP_FILE" ]; then
        # DEPS is first line of $DEP_FILE with some modifications:
        # - start with a single space
        # - the # symbol is treated as comment which is removed
        # - lower-case words "foo" are replaced by "$(inst_foo)"
        DEPS=`sed 's/^ */ /; s/ *#.*//; s/ \([a-z0-9][a-z0-9_]*\)/ $(inst_\1)/g; q' $DEP_FILE`
    elif [ "$TYPE" = optional ]; then
        DEPS=' | $(STANDARD_PACKAGES)' # default for optional packages
    elif [ "$TYPE" = script ]; then
        DEPS=' | $(STANDARD_PACKAGES)' # default for script-only packages
    elif [ "$TYPE" = pip ]; then
        DEPS=' | $(inst_pip)'
    else
        DEPS=""
    fi

    if [ "$TYPE" = pip ]; then
        # Special rules using PIP
        echo >&7 "$PKG_NAME:$DEPS"
        echo >&7 "	sage-logger -p 'sage --pip install $PKG_NAME' \$(SAGE_LOGS)/$PKG_NAME.log"
        echo >&7

        # Add a target to remove the "installed" file for "sage -f"
        echo >&7 "$PKG_NAME-clean:"
        echo >&7 "	-sage --pip uninstall -y $PKG_NAME"
        echo >&7
    elif [ "$TYPE" = script ]; then
        # Just run the spkg-install script
        echo >&7 "$PKG_NAME:$DEPS"
        echo >&7 "	cd '$SAGE_ROOT' && \\"
        echo >&7 "	    source '$SAGE_ROOT/src/bin/sage-env' && \\"
        echo >&7 "	    sage-logger -p '$SAGE_ROOT/build/pkgs/$PKG_NAME/spkg-install' \$(SAGE_LOGS)/$PKG_NAME.log"
        echo >&7

        # Just run the spkg-uninstall script if it exists
        echo >&7 "$PKG_NAME-clean:"
        echo >&7 "	-cd '$SAGE_ROOT' && \\"
        echo >&7 "	    source '$SAGE_ROOT/src/bin/sage-env' && \\"
        echo >&7 "	    '$SAGE_ROOT/build/pkgs/$PKG_NAME/spkg-uninstall'"
        echo >&7
    else
        # Normal Sage packages
        echo >&7 "\$(INST)/$PKG_VERSION:$DEPS"
        echo >&7 "	+sage-logger -p '\$(SAGE_SPKG) $PKG_VERSION' '\$(SAGE_LOGS)/$PKG_VERSION.log'"
        echo >&7

        # Add a target with just the bare package name for "sage -i"
        echo >&7 "$PKG_NAME: \$(INST)/$PKG_VERSION"
        echo >&7

        # Add a target to remove the "installed" file for "sage -f"
        echo >&7 "$PKG_NAME-clean:"
        echo >&7 "	rm -f \$(INST)/$PKG_VERSION"
        echo >&7
    fi

    PHONY="$PHONY $PKG_NAME"
done || exit 1
echo >&7 ".PHONY:$PHONY"
}

# Close the Makefile
exec 7>&-
]


dnl AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([build/make/Makefile-auto])
AC_CONFIG_MACRO_DIR([m4])

dnl Create basic directories needed for Sage
AC_CONFIG_COMMANDS(mkdirs,
    [
        SAGE_LOGS=${ac_abs_top_builddir}/logs/pkgs
        SAGE_LOCAL=${ac_abs_top_builddir}/local
        SAGE_SHARE=${ac_abs_top_builddir}/local/share
        SAGE_INST=${ac_abs_top_builddir}/local/var/lib/sage/installed
        for d in "$SAGE_LOGS" "$SAGE_LOCAL" \
            "$SAGE_LOCAL/bin" "$SAGE_LOCAL/etc" \
            "$SAGE_LOCAL/include" "$SAGE_LOCAL/lib" \
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
    ])

AC_OUTPUT()

dnl vim:syntax=m4
