#!/usr/bin/env bash
#
# Determine whether the C compiler $CC $CFLAGS supports certain
# compiler flags.  "supported" simply means here that a very basic
# program compiles.
#
# Usage: testcflags.sh [CFLAG1] [CFLAG2]...
#
# This will test all given flags one by one (together with $CFLAGS) and
# print (on stdout) a space-separated list of supported flags.  Note
# that the order is important: if CFLAG1 works, it will be added while
# testing CFLAG2.
#
# The exit code is 0 if all given flags are supported, 1 if some flag
# is not supported, 2 if some other error occurred.
#
# For example, running
#   $ testcflags.sh --foobar -Wall -march=native
# will return (on recent gcc versions)
#   -Wall -march=native
# with exit status 1 (because --foobar is not supported).
#
# A typical use would be:
#   $ CFLAGS="$CFLAGS `testcflags.sh -ffoo-bar`"
#
# It is allowed to group flags, these will then be tested together.
# For example, running
#   $ testcflags.sh '-march=native -mno-avx'
# will test these flags together.  So, the output can be either empty
# or it can be "-march=native -mno-avx", but not "-march=native".
#
#
# AUTHORS:
#
# - Jeroen Demeyer (2012-01-27): initial version (#12367)

# - Jeroen Demeyer (2012-04-11): various fixes (#12821)
#
#*****************************************************************************
#       Copyright (C) 2012 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

usage()
{
cat >&2 <<EOF
Usage: $0 [CFLAG1] [CFLAG2]...

Determine whether the C compiler \$CC \$CFLAGS supports certain
compiler flags.  "supported" simply means here that a very basic
program compiles.

This will test all given flags one by one (together with \$CFLAGS) and
print (on stdout) a space-separated list of supported flags.  Note
that the order is important: if CFLAG1 works, it will be added while
testing CFLAG2.

The exit code is 0 if all given flags are supported, 1 if some flag
is not supported, 2 if some other error occurred.
EOF
}

if [ -z "$CC" ]; then
    echo >&2 "$0: set \$CC before running this script."
    echo
    usage
    exit 2
fi

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

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
outfile=sage-testcflags-$$
cfile=$outfile.c

cat >$cfile <<EOF
/* The following function yields assembler errors on OS X 10.7 with
 * certain versions of XCode when compiling code for a machine with
 * support for AVX (Advanced Vector Extensions). */
volatile double d;
unsigned long double_to_ulong()
{
    return (unsigned long)d;
}

int main(int argc, char **argv)
{
    return 0;
}
EOF
if [ $? -ne 0 ]; then
    echo >&2 "Error: cannot write to file `pwd`/$cfile"
    exit 2
fi

status=0
NEWCFLAGS=""
for testflag in "$@"; do
    if $CC $CFLAGS $NEWCFLAGS $testflag $cfile -o $outfile 2>/dev/null; then
        # Success
        NEWCFLAGS="$NEWCFLAGS $testflag"
    else
        # Failure
        status=1
    fi
done

# Some OS X systems create a directory with debug info,
# so we really need -r here (#13945).
rm -rf $outfile*

echo $NEWCFLAGS
exit $status
