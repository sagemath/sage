#!/usr/bin/env bash
# Invoke make with -j N when called from a make running a jobserver.
OPTIONS=
case "$MAKEFLAGS" in
    *jobserver*)
        if [ -n "$SAGE_NUM_THREADS_PARALLEL" ]; then
            # Maximum number of job slots requested from the token server.
            OPTIONS="-j$SAGE_NUM_THREADS_PARALLEL"
        fi
        ;;
esac
exec $MAKE $OPTIONS "$@"
