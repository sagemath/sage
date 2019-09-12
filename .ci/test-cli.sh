#!/bin/sh

# This script gets called from CI to run minimal tests on the sagemath image.

# Usage: ./test-cli.sh IMAGE-NAME

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

echo "Checking that Sage starts and can calculate 1+1…"
# Calculate 1+1 (remove startup messages and leading & trailing whitespace)
TWO=`docker run "$1" sage -c "'print(1+1)'" | tail -1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'`
[ "x$TWO" = "x2" ]

echo "Checking that some binaries that should be distributed with Sage are on the PATH…"
# We could also run minimal tests on these but we don't yet.
# Check that Singular and GAP are present
docker run "$1" which Singular
docker run "$1" which gap
# Check that jupyter is present (for binder)
docker run "$1" which jupyter
