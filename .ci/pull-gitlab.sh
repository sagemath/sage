#!/bin/sh

# This script gets called from CI to pull the Sage docker images that were
# built during the "build" phase to pull all the connected docker daemon
# (likely a docker-in-docker.)
# This script expects a single parameter, the base name of the docker image
# such as sagemath or sagemath-dev.
# The variable $DOCKER_IMAGE is set to the full name of the pulled image;
# source this script to use it.

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

# Pull the built images from the gitlab registry and give them the original
# names they had after built.
# Note that "set -x" prints the $CI_BUILD_TOKEN here but GitLab removes it
# automatically from the log output.
docker login -u gitlab-ci-token -p $CI_BUILD_TOKEN $CI_REGISTRY
docker pull $CI_REGISTRY_IMAGE/$1:$DOCKER_TAG
export DOCKER_IMAGE="${DOCKER_NAMESPACE:-sagemath}/$1:$DOCKER_TAG"
docker tag $CI_REGISTRY_IMAGE/$1:$DOCKER_TAG $DOCKER_IMAGE
