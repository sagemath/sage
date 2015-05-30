dnl Very loosely based on configure.ac in prereq-0.3 written by William Stein.
dnl Version 0.7 written by David Kirkby, released under the GPL version 2.
dnl in January 2010

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

dnl Older versions give "undefined macro: _m4_divert_diversion", see
dnl http://trac.sagemath.org/ticket/15606#comment:19
AC_PREREQ([2.64])

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

AX_CHECK_ROOT([AC_MSG_ERROR([You cannot build Sage as root, switch to a unpriviledged user])], [])

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


dnl Check compiler versions
buggy_gcc_version1="4.0.0"
minimum_gcc_version_for_no_hassle="4.0.1"
minimum_gcc_version_for_debugging_purposes="3.4.0"

#---------------------------------------------------------
#   Process options and environment variables

# Check --enable-compiler-checks (yes by default)
AC_ARG_ENABLE([compiler-checks], [Check versions and presence of C, C++ and Fortran compilers (default: yes)],
    [enable_compiler_checks=$enableval], [enable_compiler_checks=yes])

# Import environment variables.
source src/bin/sage-env || AC_MSG_ERROR([failed to source sage-env])

#---------------------------------------------------------
# Check some programs needed actually exist.
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

# To build Python on multi-arch Debian-based systems, we need
# dpkg-architecture. Since we need dpkg-architecture to determine
# whether we're on a multi-arch system and require dpkg-architecture,
# we simply require it always on Debian-based systems.
AC_CHECK_PROG(found_dpkg, dpkg, yes, no)
AC_CHECK_PROG(found_dpkg_arch, dpkg-architecture, yes, no)
if test x$found_dpkg = xyes && test x$found_dpkg_arch = xno
then
    AC_MSG_NOTICE([You do not have 'dpkg-architecture', which is required to build])
    AC_MSG_NOTICE([Python on multi-arch Debian-based systems. This includes all recent])
    AC_MSG_NOTICE([Debian and Ubuntu systems. You can install this with:])
    AC_MSG_NOTICE([  sudo apt-get install dpkg-dev])
    AC_MSG_ERROR([Exiting, since AC_PACKAGE_NAME requires dpkg-architecture on Debian])
fi

#---------------------------------------------------------
#       C/C++/Fortran compilers

# First check for programs we need.
AC_LANG(C)
AC_PROG_CC()
AC_PROG_CPP()

AC_LANG(C++)
AC_PROG_CXX()

AC_LANG(Fortran)
AC_PROG_FC()

if test "x$CXX" = x
then
    AC_MSG_ERROR([Exiting, since a C++ compiler was not found.])
fi

if test $enable_compiler_checks = yes
then
    if test "x$CC" = x
    then
        AC_MSG_ERROR([Exiting, since a C compiler was not found.])
    fi
    if test "x$FC" = x
    then
        AC_MSG_ERROR([Exiting, since a Fortran compiler was not found.])
    fi
fi

# Next one should check for header files.
# complex.h is one that might not exist on older systems.
AC_LANG(C++)
AC_CHECK_HEADER([complex.h],[],[
    AC_MSG_ERROR([Exiting, since you do not have the 'complex.h' header file.])
])

# Next one should check for types.
# None needed

# Next one should check for structures.
# None needed

# Next one should check for compiler characterists.

# Check that we can compile C99 code
AC_LANG(C)
AC_PROG_CC_C99()
if test $enable_compiler_checks = yes
then
    if test "x$ac_cv_prog_cc_c99" = xno
    then
        AC_MSG_ERROR([Exiting, as your C compiler cannot compile C99 code])
    fi
fi

# Check that the Fortran compiler accepts free-format source code
# (as opposed to the older fixed-format style from Fortran 77).
# This helps verify the compiler works too, so if some idiot
# sets FC to /usr/bin/ls, we will at least know it's
# not a working Fortran compiler.
AC_LANG(Fortran)
if test $enable_compiler_checks = yes
then
   # see http://www.gnu.org/software/hello/manual/autoconf/Fortran-Compiler.html
   AC_FC_FREEFORM([],
   [
   AC_MSG_NOTICE([Your Fortran compiler does not accept free-format source code])
   AC_MSG_NOTICE([which means the compiler is either seriously broken, or])
   AC_MSG_NOTICE([is too old to build Sage.])
   AC_MSG_ERROR([Exiting, as the Fortran compiler is not suitable])
   ])
fi

if test $enable_compiler_checks = yes
then
    # Check that all compilers (C, C++, Fortan) are either all GNU
    # compiler or all non-GNU compilers.  If not, there is a problem, as
    # mixing GNU and non-GNU compilers is likely to cause problems.
    if test x$GCC = xyes && test x$GXX != xyes
    then
        AC_MSG_NOTICE([You are trying to use gcc but not g++])
        AC_MSG_ERROR([The mixing of GNU and non-GNU compilers is not permitted])
    fi
    if test x$GXX = xyes && test x$ac_cv_fc_compiler_gnu != xyes
    then
        AC_MSG_NOTICE([You are trying to use g++ but not gfortran])
        AC_MSG_ERROR([The mixing of GNU and non-GNU compilers is not permitted])
    fi
    if test x$ac_cv_fc_compiler_gnu = xyes && test x$GCC != xyes
    then
        AC_MSG_NOTICE([You are trying to use gfortran but not gcc])
        AC_MSG_ERROR([The mixing of GNU and non-GNU compilers is not permitted])
    fi
fi

# The following tests check the version of the compilers (if GNU)
# are all the same. If non-GNU compilers are used, then no such
# checks are performed.
if test $enable_compiler_checks = yes
then
if test x$GCC = xyes
then
    # Thank you to Andrew W. Nosenko andrew.w.nosenko@gmail.com
    # who answered my query about testing of gcc versions on
    # the autoconf@gnu.org mailing list.
    # AS_VERSION_COMPARE(ver-1, ver-2, [action-if-less], [action-if-eq], [action-if-greater])
    AX_GCC_VERSION
    AX_GXX_VERSION
    AS_VERSION_COMPARE([$GCC_VERSION], [$GXX_VERSION],
    [
    AC_MSG_NOTICE([gcc ($GCC_VERSION) and g++ ($GXX_VERSION) are not the same version])
    AC_MSG_NOTICE([which they must be. Check your setting of CC and CXX])
    AC_MSG_ERROR([Exiting, since the C and C++ compilers have different versions])
    ],[],[
    AC_MSG_NOTICE([gcc ($GCC_VERSION) and g++ ($GXX_VERSION) are not the same version])
    AC_MSG_NOTICE([which they must be. Check your setting of CC and CXX])
    AC_MSG_ERROR([Exiting, since the C and C++ compilers have different versions])])

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
        AC_MSG_NOTICE([Although gcc and g++ are both version $GCC_VERSION])
        AC_MSG_NOTICE([the Fortran compiler $FC is some other version.])
        AC_MSG_NOTICE([The output from $FC --version is below.])
        echo ""
        $FC --version 2>&1
        echo ""
        AC_MSG_ERROR([Exiting, since the Fortran compiler is not the same version as the C and C++ compilers])
    else
        AC_MSG_NOTICE([Excellent, the C, C++ and Fortran compilers are all GCC $GCC_VERSION])
    fi

    # Exit if the version of GCC is known to be too old, and so old
    # we have no intension whatsoever of trying to make Sage work with it.
    AS_VERSION_COMPARE([$GCC_VERSION], [$minimum_gcc_version_for_debugging_purposes],[
        AC_MSG_NOTICE([GCC $GCC_VERSION is too old and can not build AC_PACKAGE_NAME. ])
        AC_MSG_NOTICE([Please use a GCC of at least $minimum_gcc_version_for_no_hassle ])
        AC_MSG_NOTICE([There are no plans whatsoever to support GCC $GCC_VERSION])
        AC_MSG_ERROR([Exiting, due to the use of a version of GCC that is too old])
    ], [],[])

    # Exit if Sage is *precisely* version $buggy_gcc_version1, as that is very buggy.
    # If any later versions of GCC are found to be buggy (which would never surprise me)
    # We can easily add a buggy_gcc_version2 and repeat the next 5 lines.
    # At the time of writing (28th September 2009) that is version 4.0.0,
    # but rather than hard-code that, it is set as a variable.
    AS_VERSION_COMPARE([$GCC_VERSION], [$buggy_gcc_version1],[],[
        AC_MSG_NOTICE([GCC $buggy_gcc_version1  is very buggy and can not build AC_PACKAGE_NAME.])
        AC_MSG_NOTICE([Please use a gcc of at least $minimum_gcc_version_for_no_hassle.])
        AC_MSG_ERROR([Exiting, due to the use of a version of GCC that is too buggy])
     ],[])

    # Issue a warning if gcc is too old to work with all packages, but
    # someone wants to try getting one or more packages work with
    # an earlier gcc. At the time of writing, (28th Sept 2009), ratpoints
    # has such an issue, requiring version 4.0.1, but we would like to
    # get it to work with version 3.4.x
    AS_VERSION_COMPARE([$GCC_VERSION], [$minimum_gcc_version_for_no_hassle],[
        AC_MSG_NOTICE([******************************************************])
        AC_MSG_NOTICE([******************************************************])
        AC_MSG_NOTICE([******************************************************])
        AC_MSG_NOTICE([GCC $GCC_VERSION is too old and can not build AC_PACKAGE_NAME.])
        AC_MSG_NOTICE([Please use a gcc of at least $minimum_gcc_version_for_no_hassle])
        AC_MSG_NOTICE([if you just want AC_PACKAGE_NAME to build without problems.])
        AC_MSG_NOTICE([])
        if test "${SAGE_USE_OLD_GCC+set}" = set; then
           AC_MSG_NOTICE([Since the variable SAGE_USE_OLD_GCC was set, the])
           AC_MSG_NOTICE([build will continue, but it will fail without changes])
           AC_MSG_NOTICE([to the Sage source code. You can be 100% sure of that.])
        else
           AC_MSG_NOTICE([If you want to try building Sage with a GCC 3.4.x     ])
           AC_MSG_NOTICE([with a view to debugging the problems which stop it   ])
           AC_MSG_NOTICE([working on a gcc 3.4 series compiler, set the         ])
           AC_MSG_NOTICE([environment variable SAGE_USE_OLD_GCC to something non])
           AC_MSG_NOTICE([empty. But if you just want Sage to work, upgrade GCC ])
           AC_MSG_NOTICE([******************************************************])
           AC_MSG_NOTICE([******************************************************])
           AC_MSG_ERROR([Exiting, due to the use of a version of GCC that is too old])
        fi
    ],
[
AC_MSG_NOTICE([Good, gcc and g++ are are just new enough as GCC $minimum_gcc_version_for_no_hassle])
AC_MSG_NOTICE([is the minimum version needed needed to build AC_PACKAGE_NAME])
]
,
[
AC_MSG_NOTICE([Excellent, GCC $GCC_VERSION is later than the minimum])
AC_MSG_NOTICE([needed to build AC_PACKAGE_NAME, which is GCC version $minimum_gcc_version_for_no_hassle])
])

    # AS_VERSION_COMPARE(ver-1, ver-2, [action-if-less], [action-if-eq], [action-if-greater])
else
    AC_MSG_WARN([You have a non-GNU compiler, but AC_PACKAGE_NAME has never been built])
    AC_MSG_WARN([successfully with any compiler other than GCC, despite])
    AC_MSG_WARN([some attempts made on Solaris to use Sun's compiler, which])
    AC_MSG_WARN([produces faster code than GCC. However, the AC_PACKAGE_NAME developers])
    AC_MSG_WARN([want AC_PACKAGE_NAME to work with other compilers, so please try.])
    AC_MSG_WARN([The AC_PACKAGE_NAME developers would welcome any feedback you can give.])
    AC_MSG_WARN([Please visit http://groups.google.com/group/sage-devel])
    AC_MSG_WARN([If you just want to use AC_PACKAGE_NAME, we suggest the use of])
    AC_MSG_WARN([GCC of at least version $minimum_gcc_version_for_no_hassle])
fi
fi  # test $enable_compiler_checks = yes

# Testing for library functions
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


dnl AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([build/Makefile-auto])
AC_CONFIG_MACRO_DIR([m4])

AC_OUTPUT()

dnl vim:syntax=m4
