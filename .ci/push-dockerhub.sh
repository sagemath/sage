#!/bin/sh

# This script gets called from CI to push our docker images to
# $DOCKER_NAMESPACE/sagemath* on the Docker Hub.
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

[ -z "$DOCKER_TAG" ] && (echo "Can not push untagged build."; exit 0)

# Push the built images to the docker hub (and fail silently if
# DOCKER_USER/SECRET_DOCKER_PASS have not been configured.)
if [ -z "$DOCKER_USER" -o -z "$SECRET_DOCKER_PASS" ]; then
    echo "DOCKER_USER/SECRET_DOCKER_PASS variables have not been configured in your Continuous Integration setup. Not pushing built images to Docker Hub."
else
  cat "$SECRET_DOCKER_PASS" | docker login -u $DOCKER_USER --password-stdin
  docker push ${DOCKER_NAMESPACE:-sagemath}/$1:$DOCKER_TAG

  # For historical reasons, we also provide a -py3 tag. It's identical to the non-py3 tag.
  docker tag ${DOCKER_NAMESPACE:-sagemath}/$1:$DOCKER_TAG ${DOCKER_NAMESPACE:-sagemath}/$1:$DOCKER_TAG-py3
  docker push ${DOCKER_NAMESPACE:-sagemath}/$1:$DOCKER_TAG-py3
fi
