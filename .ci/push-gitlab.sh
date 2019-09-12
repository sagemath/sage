#!/bin/sh

# This script gets called from CI to push our docker images to registry
# configured in GitLab. (Mostly, so we can pull them again to push them to the
# Docker Hub.)
# This script expects a single parameter, the base name of the docker image
# such as sagemath or sagemath-dev.

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

# Note that "set -x" prints the $CI_BUILD_TOKEN here but GitLab removes it
# automatically from the log output.
docker login -u gitlab-ci-token -p $CI_BUILD_TOKEN $CI_REGISTRY
docker tag ${DOCKER_NAMESPACE:-sagemath}/$1:$DOCKER_TAG $CI_REGISTRY_IMAGE/$1:$DOCKER_TAG
docker push $CI_REGISTRY_IMAGE/$1:$DOCKER_TAG
