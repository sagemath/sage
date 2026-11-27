#!/bin/sh

# This script gets called from CI to establish the name of the current docker
# tag to build and also the image which is used to seed the cache.

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

# The maintainer of the CI environment, e.g., the people administrating the
# SageMath account on gitlab.com, can decide to inject an arbitrary
# base64-encoded script early into the CI execution. The script has access to
# all the CI variables, e.g., to find out which job is being executed, and
# it could do things such as "exit 1" to fail early, or "git merge" a fix for a
# known bug in CI. The CI_MONKEY_PATCH could of course also curl a more
# complicated script and execute that.
if [ -n "$CI_MONKEY_PATCH" ]; then
    SCRIPT=$(echo "$CI_MONKEY_PATCH" | base64 -d)
    $SCRIPT
fi

# From the docker documentation: "A tag name must be valid ASCII and may
# contain lowercase and uppercase letters, digits, underscores, periods and
# dashes. A tag name may not start with a period or a dash and may contain a
# maximum of 128 characters."
export DOCKER_TAG=`echo $DOCKER_TAG | tr -d '[:space:]' | tr -c '[:alnum:]_.-' '-' | sed 's/^[-.]*//' | cut -c1-128`

[[ -z "$DOCKER_TAG" ]] && export DOCKER_TAG=none
[[ "$DOCKER_TAG" = "master" ]] && export DOCKER_TAG=latest

export DOCKER_IMAGE_CLI=${DOCKER_NAMESPACE:-sagemath}/sagemath:$DOCKER_TAG
export DOCKER_IMAGE_DEV=${DOCKER_NAMESPACE:-sagemath}/sagemath-dev:$DOCKER_TAG

export DOCKER_IMAGE_BINDER="${DOCKER_NAMESPACE:-sagemath}/sagemath:${CI_COMMIT_SHA}"

# Seed the build cache with this image (set to source-clean to build from
# scratch.)
export ARTIFACT_BASE=${ARTIFACT_BASE:-$DEFAULT_ARTIFACT_BASE}
