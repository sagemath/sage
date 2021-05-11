#!/bin/sh

# Source this script during a CI run to set environment variables and print
# some informational messages about the system we are running on.

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

# CircleCI has no mechanism to hide secret variables.
# Therefore we roll our own to protect $SECRET_* variables.
. .ci/protect-secrets.sh
# Collect debug infos about the system we are running on
.ci/describe-system.sh
# Set MAKEFLAGS and SAGE_NUM_THREADS
. .ci/setup-make-parallelity.sh

# Set DOCKER_TAG according to the current branch/tag
export DOCKER_TAG=${CIRCLE_TAG:-$CIRCLE_BRANCH}
. .ci/update-env.sh
