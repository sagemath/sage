#!/bin/sh
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


# Resolve all symbolic links in a filename.  This more or less behaves
# like "readlink -f" except that it does not convert the filename to an
# absolute path (a relative path remains relative), nor does it treat
# "." or ".." specially.
#
# WARNING: this function is copy/pasted from src/bin/sage-env, and
# would deserve to be factored out if we could -- but we need it here so
# we can actually find other files in SAGE_ROOT!
# The copies should be kept synchronized.
resolvelinks() {
    # $in is what still needs to be converted (normally has no starting slash)
    in="$1"
    # $out is the part which is converted (normally ends with trailing slash)
    out="./"

    # Move stuff from $in to $out
    while [ -n "$in" ]; do
        # Normalize $in by replacing consecutive slashes by one slash
        in=$(echo "${in}" | sed 's://*:/:g')

        # If $in starts with a slash, remove it and set $out to the root
        in_without_slash=${in#/}
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
            in_without_slash=${in#/}
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
        in="${in#${f}}"
        in="${f_resolved}${in}"
    done

    # Return $out
    echo "$out"
}

# Determine SAGE_ROOT from $0. This is so we can find the src/bin/sage script.
#
# This uses "resolvelinks" to support the use case of symlinking this script
# into a directory that is in PATH, which was longstanding recommended practice.
#
# (The updated README.md and installation manual from Issue #33787 recommend to
#  symlink the installed version of the src/bin/sage script instead.)

# Get the path to $0 (this shell script) with all symbolic links resolved
SAGE_ROOT=`resolvelinks "$0"` || \
    SAGE_ROOT="$0"

# Get the directory component
SAGE_ROOT="${SAGE_ROOT%/*}"

# Make SAGE_ROOT absolute
SAGE_ROOT=`cd "$SAGE_ROOT" && pwd -P`
if [ $? -ne 0 ]; then
    echo >&2 "$0: cannot determine SAGE_ROOT directory"
    exit 1
fi

export SAGE_ROOT

# If this is a freshly-unpacked binary tarball then run the installer
# Note: relocate-once.py deletes itself upon successful completion
if [ -x "$SAGE_ROOT/relocate-once.py" ]; then
    "$SAGE_ROOT/relocate-once.py"
    if [ $? -ne 0 ]; then
        echo >&2 "Error running the script 'relocate-once.py'."
        exit 1
    fi
fi

# Run the actual Sage script
if [ -x "$SAGE_ROOT/src/bin/sage" ]; then
    exec "$SAGE_ROOT/src/bin/sage" "$@"
else
    echo >&2 "$0: no Sage installation found in \$SAGE_ROOT=$SAGE_ROOT"
    exit 1
fi
