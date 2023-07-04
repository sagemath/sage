dnl A function to test for the existence and usability of particular
dnl versions of the PPL, defining macros containing the required paths.
dnl Copyright (C) 1997 Owen Taylor
dnl Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
dnl Copyright (C) 2010-2016 BUGSENG srl (http://bugseng.com)
dnl
dnl This file is part of the Parma Polyhedra Library (PPL).
dnl
dnl The PPL is free software; you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by the
dnl Free Software Foundation; either version 3 of the License, or (at your
dnl option) any later version.
dnl
dnl The PPL is distributed in the hope that it will be useful, but WITHOUT
dnl ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
dnl FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software Foundation,
dnl Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02111-1307, USA.
dnl
dnl For the most up-to-date information see the Parma Polyhedra Library
dnl site: http://bugseng.com/products/ppl/ .

dnl AM_PATH_PPL([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for PPL, and define PPL_CPPFLAGS, PPL_LDFLAGS, ... what else?

AC_DEFUN([AM_PATH_PPL],
[
dnl Get the required information from the ppl-config program.
AC_ARG_WITH(ppl-prefix,
  AS_HELP_STRING([--with-ppl-prefix=PREFIX],
    [prefix used to configure the PPL]),
  ppl_prefix="$withval",
  ppl_prefix="")
AC_ARG_WITH(ppl-exec-prefix,
  AS_HELP_STRING([--with-ppl-exec-prefix=PREFIX],
    [exec-prefix used to configure the PPL]),
  ppl_exec_prefix="$withval",
  ppl_exec_prefix="")
AC_ARG_ENABLE(ppl-test,
  AS_HELP_STRING([--disable-ppltest],
    [do not try to compile and run a test PPL program]),
  ,
  enable_ppltest=yes)

if test "x$ppl_exec_prefix" != x
then
  ppl_config_args="$ppl_config_args --exec-prefix=$ppl_exec_prefix"
  if test "x${PPL_CONFIG+set}" != xset
  then
    PPL_CONFIG="$ppl_exec_prefix/bin/ppl-config"
  fi
fi
if test "x$ppl_prefix" != x
then
  ppl_config_args="$ppl_config_args --prefix=$ppl_prefix"
  if test "x${PPL_CONFIG+set}" != xset
  then
    PPL_CONFIG="$ppl_prefix/bin/ppl-config"
  fi
fi

AC_PATH_PROG(PPL_CONFIG, ppl-config, no)
min_ppl_version=ifelse([$1], ,0.0,$1)
if test \( "x$min_ppl_version" = "x0.0" \) -o \( "x$min_ppl_version" = "x0.0.0" \)
then
  AC_MSG_CHECKING([for the Parma Polyhedra Library])
else
  AC_MSG_CHECKING([for the Parma Polyhedra Library, version >= $min_ppl_version])
fi
no_ppl=""
if test $PPL_CONFIG = no
then
  no_ppl=yes
else
  PPL_CPPFLAGS=`$PPL_CONFIG $ppl_config_args --cppflags`
  PPL_LDFLAGS=""
  PPL_LIBS=""
  for flag in $($PPL_CONFIG $ppl_config_args --ldflags); do
    dnl Check if each "ldflag" starts with -l or not. The ones that do
    dnl start with -l are libs and belong in PPL_LIBS, the others belong
    dnl in PPL_LDFLAGS. The LDFLAGS and LIBS variables get appended at
    dnl different locations in the link command, so the distinction is
    dnl not academic.
    if test "x${flag#-l}" = "x$flag"; then
      dnl this flag doesn't start with -l
      PPL_LDFLAGS="$PPL_LDFLAGS $flag"
    else
      PPL_LIBS="$PPL_LIBS $flag"
    fi
  done
  ppl_config_version="`$PPL_CONFIG $ppl_config_args --version`"

  if test "x$enable_ppltest" = xyes
  then
    ac_save_CPPFLAGS="$CPPFLAGS"
    ac_save_LIBS="$LIBS"
    CPPFLAGS="$CPPFLAGS $PPL_CPPFLAGS"
    LIBS="$LIBS $PPL_LIBS"

dnl Now check if the installed PPL is sufficiently new.
dnl (Also sanity checks the results of ppl-config to some extent.)

    AC_LANG_PUSH(C++)

    rm -f conf.ppltest
    AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <ppl.hh>
#include <iostream>
#include <cstdio>
#include <cstdlib>

namespace PPL = Parma_Polyhedra_Library;

using std::cout;
using std::endl;

int
main() {
  system("touch conf.ppltest");

  unsigned min_ppl_major, min_ppl_minor, min_ppl_revision, min_ppl_beta;
  int n = sscanf("$min_ppl_version",
                 "%u.%u.%upre%u%*c",
                 &min_ppl_major, &min_ppl_minor,
                 &min_ppl_revision, &min_ppl_beta);
  bool min_ppl_version_ok = true;
  if (n == 4) {
    if (min_ppl_beta == 0)
      min_ppl_version_ok = false;
  }
  else if (n == 3) {
    n = sscanf("$min_ppl_version",
               "%u.%u.%u%*c",
               &min_ppl_major, &min_ppl_minor, &min_ppl_revision);
    if (n != 3)
      min_ppl_version_ok = false;
    else
      min_ppl_beta = 0;
  }
  else if (n == 2) {
    n = sscanf("$min_ppl_version",
               "%u.%upre%u%*c",
               &min_ppl_major, &min_ppl_minor, &min_ppl_beta);
    if (n == 3) {
      if (min_ppl_beta == 0)
        min_ppl_version_ok = false;
      else
        min_ppl_revision = 0;
    }
    else if (n == 2) {
      n = sscanf("$min_ppl_version",
                 "%u.%u%*c",
                 &min_ppl_major, &min_ppl_minor);
      if (n != 2)
        min_ppl_version_ok = false;
      else {
        min_ppl_revision = 0;
        min_ppl_beta = 0;
      }
    }
    else
      min_ppl_version_ok = false;
  }
  else
    min_ppl_version_ok = false;

  if (!min_ppl_version_ok) {
    cout << "illegal version string '$min_ppl_version'"
         << endl;
    return 1;
  }

  if (strcmp("$ppl_config_version", PPL::version()) != 0) {
    cout << "\n*** 'ppl-config --version' returned $ppl_config_version, "
            "but PPL version "
         << PPL::version()
         << "\n*** was found!  If ppl-config was correct, then it is best"
            "\n*** to remove the old version of PPL."
            "  You may also be able to fix the error"
            "\n*** by modifying your LD_LIBRARY_PATH environment variable,"
            " or by editing"
            "\n*** /etc/ld.so.conf."
            "  Make sure you have run ldconfig if that is"
            "\n*** required on your system."
            "\n*** If ppl-config was wrong, set the environment variable"
            " PPL_CONFIG"
            "\n*** to point to the correct copy of ppl-config,"
            " and remove the file config.cache"
            "\n*** before re-running configure."
         << endl;
      return 1;
  }
  else if (strcmp(PPL_VERSION, PPL::version()) != 0) {
    cout << "\n*** PPL header file (version " PPL_VERSION ") does not match"
         << "\n*** library (version " << PPL::version() << ")"
         << endl;
      return 1;
  }
  else if (PPL_VERSION_MAJOR < min_ppl_major
           || (PPL_VERSION_MAJOR == min_ppl_major
              && PPL_VERSION_MINOR < min_ppl_minor)
           || (PPL_VERSION_MAJOR == min_ppl_major
              && PPL_VERSION_MINOR == min_ppl_minor
              && PPL_VERSION_REVISION < min_ppl_revision)
           || (PPL_VERSION_MAJOR == min_ppl_major
              && PPL_VERSION_MINOR == min_ppl_minor
              && PPL_VERSION_REVISION == min_ppl_revision
              && PPL_VERSION_BETA < min_ppl_beta)) {
      cout << "\n*** An old version of PPL (" PPL_VERSION ") was found."
              "\n*** You need at least PPL version $min_ppl_version."
              "  The latest version of"
              "\n*** PPL is always available from ftp://ftp.cs.unipr.it/ppl/ ."
              "\n***"
              "\n*** If you have already installed a sufficiently new version,"
              " this error"
              "\n*** probably means that the wrong copy of the ppl-config"
              " program is"
              "\n*** being found.  The easiest way to fix this is to remove"
              " the old version"
              "\n*** of PPL, but you can also set the PPL_CONFIG environment"
              " variable to point"
              "\n*** to the correct copy of ppl-config.  (In this case,"
              " you will have to"
              "\n*** modify your LD_LIBRARY_PATH environment"
              " variable or edit /etc/ld.so.conf"
              "\n*** so that the correct libraries are found at run-time.)"
           << endl;
      return 1;
  }
  return 0;
}
]])],[],[no_ppl=yes],[echo $ac_n "cross compiling; assumed OK... $ac_c"])

    AC_LANG_POP

    CPPFLAGS="$ac_save_CPPFLAGS"
    LIBS="$ac_save_LIBS"
  fi
fi

if test "x$no_ppl" = x
then
  AC_MSG_RESULT(yes)
  ifelse([$2], , :, [$2])
else
  AC_MSG_RESULT(no)
  if test x"$PPL_CONFIG" = xno
  then
    echo "*** The ppl-config script installed by PPL could not be found."
    echo "*** If the PPL was installed in PREFIX, make sure PREFIX/bin is in"
    echo "*** your path, or set the PPL_CONFIG environment variable to the"
    echo "*** full path to ppl-config."
  else
    if test -f conf.ppltest
    then
      :
    else
      echo "*** Could not run PPL test program, checking why..."
      CPPFLAGS="$CPPFLAGS $PPL_CPPFLAGS"
      LIBS="$LIBS $PPL_LIBS"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <ppl.hh>
using namespace Parma_Polyhedra_Library;
]], [[
  return version_major() || version_minor()
  || version_revision() || version_beta();
]])],[
  echo "*** The test program compiled, but did not run.  This usually means"
  echo "*** that the run-time linker is not finding the PPL or finding the"
  echo "*** wrong version of the PPL.  If it is not finding the PPL, you will"
  echo "*** need to set your LD_LIBRARY_PATH environment variable, or edit"
  echo "*** /etc/ld.so.conf to point to the installed location.  Also, make"
  echo "*** sure you have run ldconfig if that is required on your system."
  echo "***"
  echo "*** If you have an old version installed, it is best to remove it,"
  echo "*** although you may also be able to get things to work by modifying"
  echo "*** LD_LIBRARY_PATH."
],[
  echo "*** The test program failed to compile or link. See the file"
  echo "*** config.log for the exact error that occured.  This usually means"
  echo "*** the PPL was incorrectly installed or that someone moved the PPL"
  echo "*** since it was installed.  In both cases you should reinstall"
  echo "*** the library."
])
      CPPFLAGS="$ac_save_CPPFLAGS"
      LIBS="$ac_save_LIBS"
    fi
  fi
  PPL_CPPFLAGS=""
  PPL_LDFLAGS=""
  PPL_LIBS=""
  ifelse([$3], , :, [$3])
fi
AC_SUBST(PPL_CPPFLAGS)
AC_SUBST(PPL_LDFLAGS)
AC_SUBST(PPL_LIBS)
rm -f conf.ppltest
])
