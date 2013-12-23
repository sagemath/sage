#!/usr/bin/env bash
#
# Sage: a free open-source mathematics software system
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


# Set SAGE_ROOT to the location of the sage install, i.e. the directory
# containing this shell script.  If unset, we will try to figure it out
# automatically.
#SAGE_ROOT=/path/to/sage-version



# Resolve all symbolic links in a filename.  This more or less behaves
# like "readlink -f" except that it does not convert the filename to an
# absolute path (a relative path remains relative), nor does it treat
# "." or ".." specially.  See Trac ticket #5852.
resolvelinks() {
    # $in is what still needs to be converted (normally has no starting slash)
    in="$1"
    # $out is the part which is converted (normally ends with trailing slash)
    out="./"

    # Move stuff from $in to $out
    while [ -n "$in" ]; do
        # Normalize $in by replacing consecutive slashes by one slash
        while { in_single_slash=${in//\/\//\/}; [ "$in" != "$in_single_slash" ]; }; do
            in=$in_single_slash
        done

        # If $in starts with a slash, remove it and set $out to the root
        in_without_slash=${in/#\//}
        if [ "$in" != "$in_without_slash" ]; then
            in=$in_without_slash
            out="/"
            continue
        fi

        # Check that the directory $out exists by trying to cd to it.
        # If this fails, then cd will show an error message (unlike
        # test -d "$out"), so no need to be more verbose.
        ( cd "$out" ) || return $?


        # Get the first component of $in
        f=${in%%/*}

        # If it is not a symbolic link, simply move it to $out
        if [ ! -L "$out$f" ]; then
            in=${in#"$f"}
            out="$out$f"

            # If the new $in starts with a slash, move it to $out
            in_without_slash=${in/#\//}
            if [ "$in" != "$in_without_slash" ]; then
                in=$in_without_slash
                out="$out/"
            fi
            continue
        fi

        # Now resolve the symbolic link "$f"
        f_resolved=`readlink -n "$out$f" 2>/dev/null`
        status=$?
        # status 127 means readlink could not be found.
        if [ $status -eq 127 ]; then
            # We don't have "readlink", try a stupid "ls" hack instead.
            # This will fail if we have filenames like "a -> b".
            fls=`ls -l "$out$f" 2>/dev/null`
            status=$?
            f_resolved=${fls##*-> }

            # If $fls equals $f_resolved, then certainly
            # something is wrong
            if [ $status -eq 0 -a "$fls" = "$f_resolved" ]; then
                echo >&2 "Cannot parse output from ls -l '$out$f'"
                return 1
            fi
        fi
        if [ $status -ne 0 ]; then
            echo >&2 "Cannot read symbolic link '$out$f'"
            return $status
        fi

        # In $in, replace $f by $f_resolved (leave $out alone)
        in=${in/#"$f"/"$f_resolved"}
    done

    # Return $out
    echo "$out"
}

# If SAGE_ROOT is not given, find it out from $0
if [ -z "$SAGE_ROOT" ];  then
    # Get the path to $0 (this shell script) with all symbolic links
    # resolved
    SAGE_ROOT=`resolvelinks "$0"` || \
    SAGE_ROOT="$0"

    # Get the directory component
    SAGE_ROOT="${SAGE_ROOT%/*}"
fi

# Make SAGE_ROOT absolute
SAGE_ROOT=`cd "$SAGE_ROOT" && pwd -P`
if [ $? -ne 0 ]; then
    echo >&2 "$0: cannot determine SAGE_ROOT directory"
    exit 1
fi

export SAGE_ROOT

# Run the actual Sage script
if [ -x "$SAGE_ROOT/src/bin/sage" ]; then
    "$SAGE_ROOT/src/bin/sage" "$@"
elif [ -x "$SAGE_ROOT/local/bin/sage" ]; then # if in a stripped binary
    "$SAGE_ROOT/local/bin/sage" "$@"
else
    echo >&2 "$0: no Sage installation found in \$SAGE_ROOT=$SAGE_ROOT"
    exit 1
fi

# Kill all processes in the current process group.  In practice, this
# means all child processes not running their own pty.
# Uncomment this if you have trouble with orphans.
# We do this in the background after waiting 10 seconds to give the
# various processes (including this shell script) some time to exit
# cleanly.
# { sleep 10; kill -9 -$$ 2>/dev/null; } &
