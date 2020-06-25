#!/bin/sh
# Determine the type of C compiler, which can later
# be used to determine the flags the C compiler
# will want. This is done by testing what the
# C compiler's pre-processor defines. For example,
# gcc will always define __GNUC__
#
# This will typically be called as (don't quote $CC):
# testcc.sh $CC

# Copyright Dr. David Kirkby
# Released under the GPL version 2, or any later version
# the user wishes to use.
# Considerably cleaned up by Peter Jeremy.
# Further clean-up by Jeroen Demeyer (#12821).
# Helpful comments by William Stein.

# Some of the compilers have been tested, though some have
# not, due to a lack of hardware or software.

# Documentation on the compilers is taken from many source
# in particular, for the commercial Unix compilers:

# HP-UX C and C++ compiler.
# http://docs.hp.com/en/7730/newhelp0610/preprocess.htm

# IBM Compiler Reference - XL C/C++ for AIX, V10.1
# http://www-01.ibm.com/support/docview.wss?uid=swg27012860&aid=1

# Using HP C++ for Tru64 UNIX and Linux Alpha
# http://h30097.www3.hp.com/cplus/ugu_impl.html#implem_chap


# Define a function to display usage information.
usage()
{
cat <<EOF 1>&2
Usage: $0 name_or_path_to_the_C_compiler [optional arguments]

Typically, this will be called as
$0 \$CC

  The script will print one of the following to indicate the C compiler

  GCC                 - For gcc or a gcc-like C compiler on any platform
  Sun_Studio          - For Sun Studio or earlier Sun C compiler
  HP_on_Tru64         - For a C compiler produced by HP/Compaq for Tru64
  HP_on_HP-UX         - For a C compiler produced by HP/Compaq for HP-UX
  IBM_on_AIX          - For a C compiler produced by IBM for AIX
  HP_on_Alpha_Linux   - For a C compiler produced by HP for Alpha linux
                        (This script has not been tested on Alpha Linux,
                        but is based on HP's documentation)
  Unknown             - If the C compiler type is unknown

EOF
}

# Exit if the user does not supply any command line arguments.
if [ $# -lt 1 ] ; then
    usage
    exit 2
fi

# Create a test file. It does not need to be a complete
# C file, as it is only pre-processed. So there is no
# need for a 'main'

if [ -z "$SAGE_ROOT" ]; then
    echo "The SAGE_ROOT environment variable must be set to"
    echo "the root of the Sage installation"
    exit 1
fi

mkdir -p "$SAGE_LOCAL/var/tmp/sage/build"
if [ $? -ne 0 ]; then
    echo "Error while trying to create the build directory."
    exit 1
fi

cd "$SAGE_LOCAL/var/tmp/sage/build"
if [ $? -ne 0 ]; then
    echo "Error while trying to change into the build directory."
    exit 1
fi
TESTFILE=sage-testcc-$$.c

cat >$TESTFILE <<"E*O*F"
#ifdef __GNUC__
GCC
#define KNOWN_CC
#endif

#ifdef __SUNPRO_C
Sun_Studio
#define KNOWN_CC
#endif

/*
 * I've not found any official documentation to suggest __DECC
 * would be defined, but __DECCC is defined for C++,
 * and a Google on __DECC shows many hits. I do not have
 * easy access to a Tru64 system.
 */
#ifdef __digital__
#ifdef __DECC
HP_on_Tru64
#define KNOWN_CC
#endif
#endif

/* Untested, as I have no access to a Alpha Linux system. */
#ifdef __linux__
#ifdef __DECCC
HP_on_Alpha_Linux
#define KNOWN_CC
#endif
#endif

#ifdef __HP_cc
HP_on_HP-UX
#define KNOWN_CC
#endif

#ifdef __xlC__
IBM_on_AIX
#define KNOWN_CC
#endif

#ifndef KNOWN_CC
Unknown
#endif
E*O*F

"$@" -E $TESTFILE | sed -n '/^[A-Z]/{s/ //g;p;}'
rm -f $TESTFILE
exit 0
