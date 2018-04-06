#!/bin/dash

# Author:
# * Jean-Pierre Flori <jean-pierre.flori@ssi.gouv.fr>
# * Gary Zablackis <gzabl@yahoo.com>
# * Dmitrii Pasechnik <dimpase@gmail.com>
# * Erik M. Bray <erik.m.bray@gmail.com>
#
# Rebase all dlls in the SAGE_LOCAL directory (and its subdirectories),
# but do not touch the ones already stored in the system database,
# and do not update it.
# Note that subsequent calls to 'rebaseall' will not update the Sage dlls.
#
# Usage:
#
# sage-rebase.sh [--all] [sage_local] [-- additional_flags]
#
# Positional arguments:
#
# sage_local  optionally, provide the path to the $SAGE_LOCAL directory to
#             search for DLLs to rebase; otherwise the current working
#             directory is assumed to be $SAGE_LOCAL unless $SAGE_LOCAL
#             is already set in the environment
#
# Optional arguments:
#
# --all       run rebaseall instead of rebase (originally the call to
#             rebaseall was in the sage-rebaseall.sh script, but now that is
#             just a wrapper around this script)
#
# --          additional arguments passed in after -- are passed to the
#             rebase/rebaseall call in addition to the default arguments
#             passed in by this script
#
# Invoke this script from a shell after going to the SAGE_LOCAL directory.
ALL=0
REBASEFLAGS=""

while [ $# -gt 0 ]; do
    case "$1" in
        --all)
            ALL=1
            ;;
        --)
            shift
            REBASEFLAGS="$REBASEFLAGS $1"
            ;;
        *)
            if [ -z "$REBASEFLAGS" ]; then
                SAGE_LOCAL="${1%/}"
            else
                REBASEFLAGS="$REBASEFLAGS $1"
            fi
            ;;
    esac
    shift
done

if [ -z "$SAGE_LOCAL" ]; then
    # Assume we are in $SAGE_LOCAL by default (the old behavior of this script)
    SAGE_LOCAL=.
fi

FINDFLAGS="-type f ( -name *.dll -o -name *.so -o -name *.fas ) -print"
FINDFLAGS="$FINDFLAGS -o -path "$SAGE_LOCAL"/var/tmp -prune"

echo "Getting list of dlls. This may take a while..."
/bin/find "$SAGE_LOCAL" $FINDFLAGS > /tmp/sage-dlls.lst
echo "Now rebasing..."

if [ $ALL -eq 0 ]; then
    /bin/rebase -O -T /tmp/sage-dlls.lst $REBASEFLAGS
else
    /bin/rebaseall -s dll -s exe -s so -s fas -T /tmp/sage-dlls.lst $REBASEFLAGS
fi
