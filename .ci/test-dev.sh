#!/bin/sh

# This script gets called from CI to run minimal tests on the sagemath-dev image.
# This script expects a single argument, the full name of the docker image to
# test.

# Usage: ./test-dev.sh IMAGE-NAME

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set -ex

IMAGE="$1"

. .ci/setup-make-parallelity.sh

# Usage: timed_run limit args
# Runs $IMAGE with args and check that it terminates with a zero exit code in at most limit seconds.
timed_run() {
    START=`date +%s`
    docker run -e MAKEFLAGS="$MAKEFLAGS_DOCBUILD" -e SAGE_NUM_THREADS="$SAGE_NUM_THREADS_DOCBUILD" "$IMAGE" "$2"
    END=`date +%s`
    TOTAL=$((END-START))
    echo "Checking whether running \"$2\" was fast…"
    [ "$TOTAL" -lt "$1" ]
}

# Most setup should be done in well under 120 seconds but some CI machines tend
# to be really slow so we better give them some space here. Better miss a
# regression at first than create a lot of noise.
timed_run 120 true # runs make build
# Building the documentation is quite slow at the moment:
# Currently, this detects some dependencies as changed that have not changed.
# The parser in Sphinx fails to parse some .py files and adds the (more
# recently modified) .pyc files as dependencies instead. (Have a look the
# changeset that introduced this comment for more details.)
# Parts of the docbuild do not seem to be run in parallel at the moment.
timed_run $(( 180 + 1200/$SAGE_NUM_THREADS_DOCBUILD )) make # runs make build and then make
