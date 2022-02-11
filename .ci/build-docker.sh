#!/bin/sh

# This script gets called from CI to build several flavours of docker images
# which contain Sage.

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

# We speed up the build process by copying built artifacts from ARTIFACT_BASE
# during docker build. See /docker/Dockerfile for more details.
ARTIFACT_BASE=${ARTIFACT_BASE:-sagemath/sagemath-dev:develop}

# Seed our cache with $ARTIFACT_BASE if it exists.
docker pull "$ARTIFACT_BASE" > /dev/null || true

docker_build() {
    # Docker's --cache-from does not really work with multi-stage builds: https://github.com/moby/moby/issues/34715
    # So we just have to rely on the local cache.
    time docker build -f docker/Dockerfile \
--build-arg "MAKEFLAGS=${MAKEFLAGS}" --build-arg "SAGE_NUM_THREADS=${SAGE_NUM_THREADS}" --build-arg "MAKEFLAGS_DOCBUILD=${MAKEFLAGS}" --build-arg "SAGE_NUM_THREADS_DOCBUILD=${SAGE_NUM_THREADS_DOCBUILD}" --build-arg ARTIFACT_BASE=$ARTIFACT_BASE $@
}

# We use a multi-stage build /docker/Dockerfile. For the caching to be
# effective, we populate the cache by building the run/build-time-dependencies
# and the make-all target. (Just building the last target is not enough as
# intermediate targets could be discarded from the cache [depending on the
# docker version] and therefore the caching would fail for our actual builds
# below.)
docker_build --target run-time-dependencies --tag run-time-dependencies:$DOCKER_TAG .
docker_build --target build-time-dependencies --tag build-time-dependencies:$DOCKER_TAG .
docker_build --target make-all --tag make-all:$DOCKER_TAG .

# Copy docs out of the docker image to save them into browseable GitLab artifacts
container=$(docker create make-all:$DOCKER_TAG)
docker cp $container:/home/sage/sage/local/share/doc/sage/html html_
# GitLab's browser does not like symlinks, so we flatten them
rsync html_/ html/ -a --copy-links

# Build the release image without build artifacts.
docker_build --target sagemath --tag "$DOCKER_IMAGE_CLI" .
# Tag the sagemath:$DOCKER_TAG image that CI has just built as
# sagemath:$COMMIT_HASH so we can refer to it uniquely later.
docker tag "$DOCKER_IMAGE_CLI" "$DOCKER_IMAGE_BINDER"
# Display the layers of this image
docker history "$DOCKER_IMAGE_CLI"
# Build the developer image with the build artifacts intact.
# Note: It is important to build the dev image last because it might be tagged as ARTIFACT_BASE.
docker_build --target sagemath-dev --tag "$DOCKER_IMAGE_DEV" .
# Display the layers of this image
docker history "$DOCKER_IMAGE_DEV"
