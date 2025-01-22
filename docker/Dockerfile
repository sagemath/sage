################################################################################
# SageMath images for Docker                                                   #
################################################################################
# This is a description of the layout of this Dockerfile; for details on the   #
# created docker images, see the README.md please.                             #
#                                                                              #
# This Dockerfile builds sagemath (for end-users) and sagemath-dev (for        #
# developers.) It consists of lots of intermediate targets, mostly to shrink   #
# the resulting images but also to make this hopefully easier to maintain.     #
# The aims of this Dockerfile are:                                             #
# (1) Make it build in reasonable time.                                        #
# (2) It should be self-contained and work on its own, i.e., just by invoking  #
# docker build without any external orchestration script.                      #
#                                                                              #
# The idea to achieve (1) is to reuse the build artifacts from the latest      #
# develop build. This is slightly against the philosophy of a Dockerfile (which#
# should produce perfectly reproducible outputs) but building Sage from scratch#
# just takes too long at the moment to do this all the time. ARTIFACT_BASE     #
# controls which build artifacts are used. You probably want to set this to    #
# sagemath/sagemath-dev:develop which takes the latest build from the official #
# develop branch. The default is source-clean which builds Sage from scratch.  #
# If you want to understand how this works, have a look at source-from-context #
# which merges ARTIFACT_BASE with the context, i.e., the contents of the sage  #
# source directory.                                                            #
################################################################################

################################################################################
# HOWTO use this file for local builds                                         #
################################################################################
# If you don't mind downloading a 2GB docker image from time to time, you      #
# could use this file for local Sage development. As of early 2018 each build  #
# takes about five minutes but you don't have to go through the sadly frequent #
# rebuilds the whole Sage distribution...                                      #
# To build Sage, run this command from your sage/ directory:                   #
# $ docker build --build-arg MAKEFLAGS="-j4" --build-arg SAGE_NUM_THREADS="4" --build-arg ARTIFACT_BASE="sagemath/sagemath-dev:develop" -f docker/Dockerfile --target=make-build --tag sage .
# To run Sage:                                                                 #
# $ docker run -it sage                                                        #
# To run doctests:                                                             #
# $ docker run -e "MAKEFLAGS=-j4" -e "SAGE_NUM_THREADS=4" -it sage sage -tp src/sage
# Make sure that you always have the latest develop branch merged into your    #
# local branch for this to work.                                               #
################################################################################

################################################################################
# HOWTO use this file to manually create new images for Docker Hub             #
################################################################################
# Make sure you have the tag checked out that you are trying to build. Then    #
# from the root directory of the sage repository run:                          #
#                                                                              #
# bash  # open a subshell so 'set -e' does not close your shell on error       #
# export ARTIFACT_BASE=source-clean                                            #
# export DOCKER_TAG=9.6  # replace with the version you are building for       #
# export CI_COMMIT_SHA=`git rev-parse --short HEAD`                            #
# export CPUTHREADS=8                                                          #
# export RAMTHREADS=8                                                          #
# export RAMTHREADS_DOCBUILD=8                                                 #
# . .ci/update-env.sh                                                          #
# . .ci/setup-make-parallelity.sh                                              #
# .ci/build-docker.sh                                                          #
#                                                                              #
# Now you can push the relevant images to the Docker Hub:                      #
#                                                                              #
# docker push sagemath/sagemath:$DOCKER_TAG                                    #
# docker tag sagemath/sagemath:$DOCKER_TAG sagemath/sagemath:latest            #
# docker push sagemath/sagemath:latest                                         #
# docker push sagemath/sagemath-dev:$DOCKER_TAG                                #
# docker tag sagemath/sagemath-dev:$DOCKER_TAG sagemath/sagemath-dev:latest    #
# docker push sagemath/sagemath-dev:latest                                     #
################################################################################


ARG ARTIFACT_BASE=source-clean
ARG MAKE_BUILD=make-build

################################################################################
# Image containing the run-time dependencies for Sage                          #
################################################################################
FROM ubuntu:jammy AS run-time-dependencies
LABEL maintainer="Erik M. Bray <erik.bray@lri.fr>, Julian Rüth <julian.rueth@fsfe.org>, Sebastian Oehms <seb.oehms@gmail.com>"
# Set sane defaults for common environment variables.
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV SHELL /bin/bash
# Create symlinks for sage and sagemath - we copy a built sage to the target of these symlinks later.
ARG SAGE_ROOT=/home/sage/sage
RUN ln -s "$SAGE_ROOT/sage" /usr/bin/sage
RUN ln -s /usr/bin/sage /usr/bin/sagemath
# Sage needs the fortran libraries at run-time because we do not build gfortran
# with Sage but use the system's.
# We need gcc/g++ and libstdc++-10-dev to allow compilation of cython at run-time from the notebook.
# We also install sudo for the sage user, see below.
RUN apt-get -qq update \
    && apt-get -qq install -y --no-install-recommends gfortran gcc g++ libstdc++-10-dev sudo openssl \
    && apt-get -qq clean \
    && rm -r /var/lib/apt/lists/*
# Sage refuses to build as root, so we add a "sage" user that can sudo without a password.
# We also want this user at runtime as some commands in sage know about the user that was used during build.
ARG HOME=/home/sage
RUN adduser --quiet --shell /bin/bash --gecos "Sage user,101,," --disabled-password --home "$HOME" sage \
    && echo "sage ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/01-sage \
    && chmod 0440 /etc/sudoers.d/01-sage
# Run everything from now on as the sage user in sage's home
USER sage
ENV HOME $HOME
WORKDIR $HOME

################################################################################
# Image containing everything so that a make in a clone of the Sage repository #
# completes without errors                                                     #
################################################################################
FROM run-time-dependencies AS build-time-dependencies
# Install the build time dependencies & git & rdfind
RUN sudo apt-get -qq update \
    && sudo apt-get -qq install -y wget build-essential automake m4 dpkg-dev python3 libssl-dev git rdfind \
    && sudo apt-get -qq clean \
    && sudo rm -r /var/lib/apt/lists/*

################################################################################
# Image with an empty git repository in $SAGE_ROOT.                            #
################################################################################
FROM build-time-dependencies AS source-clean
ARG SAGE_ROOT=/home/sage/sage
RUN mkdir -p "$SAGE_ROOT"
WORKDIR $SAGE_ROOT
RUN git init
RUN git remote add upstream https://github.com/sagemath/sage.git

################################################################################
# Image with the build context added, i.e., the directory from which `docker   #
# build` has been called in a separate directory so we can copy files from     #
# there.                                                                       #
# This blows up the size of this docker image significantly, but we only use   #
# this image to create artifacts for our final image.                          #
# Warning: If you set ARTIFACT_BASE to something else than source-clean, the   #
# build is not going to use the build-time-dependencies target but rely on     #
# whatever tools are installed in ARTIFACT_BASE.                               #
################################################################################
FROM $ARTIFACT_BASE AS source-from-context
WORKDIR $HOME
COPY --chown=sage:sage . sage-context
# Checkout the commit that is checked out in $HOME/sage-context
# This is a bit complicated because our local .git/ is empty and we want to
# make sure that we only change the mtimes of a minimal number of files.
# 1) Restore the git checkout ARTIFACT_BASE was built from, recorded in
#    docker/.commit. (Or go directly to FETCH_HEAD if there is no history to
#    restore, i.e., set ARTIFACT_BASE=source-clean if you want to build from
#    scratch.)
# 2) Merge in FETCH_HEAD but only if it is a fast-forward, i.e., if it is an
#    ancestor of the commit restored in 1. If we would not do that we would get
#    a new commit hash in docker/.commit that is not known outside of this build
#    run. Since docker/.commit was in the history of FETCH_HEAD this should
#    automatically be a fast-forward.
# 3) Trash .git again to save some space.
ARG SAGE_ROOT=/home/sage/sage
WORKDIR $SAGE_ROOT
# We create a list of all files present in the artifact-base (with a timestamp
# of now) so we can find out later which files were added/changed/removed.
RUN find . \( -type f -or -type l \) > $HOME/artifact-base.manifest
RUN git fetch --update-shallow "$HOME/sage-context" HEAD \
    && if [ -e docker/.commit ]; then \
          git reset `cat docker/.commit` \
          || ( echo "Could not find commit `cat docker/.commit` in your local Git history. Please merge in the latest built develop branch to fix this: git fetch upstream && git merge `cat docker/.commit`." && exit 1 ) \
       else \
          echo "You are building from $ARTIFACT_BASE which has no docker/.commit file. That's a bug unless you are building from source-clean or something similar." \
          && git reset FETCH_HEAD \
          && git checkout -f FETCH_HEAD; \
       fi \
    && git merge --ff-only FETCH_HEAD \
    && git log -1 --format=%H > docker/.commit \
    && rm -rf .git
# Copy over all the untracked/staged/unstaged changes from sage-context. This
# is relevant for non-CI invocations of this Dockerfile.
WORKDIR $HOME/sage-context
RUN if git status --porcelain | read CHANGES; then \
        git -c user.name=docker-build -c user.email=docker-build@sage.invalid stash -u \
        && git stash show -p > "$HOME"/sage-context.patch; \
    else \
        touch "$HOME"/sage-context.patch; \
    fi
WORKDIR $SAGE_ROOT
RUN patch -p1 < "$HOME"/sage-context.patch

################################################################################
# Image with a built sage but without sage's documentation.                    #
################################################################################
FROM source-from-context AS make-build
# Make sure that the result runs on most CPUs.
ENV SAGE_FAT_BINARY yes
# Just to be sure Sage doesn't try to build its own GCC (even though
# it shouldn't with a recent GCC package from the system and with gfortran)
ENV SAGE_INSTALL_GCC no
# Set MAKEFLAGS and SAGE_NUM_THREADS to build things in parallel during the
# docker build. Note that these do not leak into the sagemath and sagemath-dev
# images.
ARG MAKEFLAGS="-j2"
ENV MAKEFLAGS $MAKEFLAGS
ARG SAGE_NUM_THREADS="2"
ENV SAGE_NUM_THREADS $SAGE_NUM_THREADS
RUN ./bootstrap
# Old default before https://github.com/sagemath/sage/issues/32406
RUN ./configure --disable-editable
RUN make build

################################################################################
# Image with a full build of sage and its documentation.                       #
################################################################################
FROM $MAKE_BUILD AS make-all
# The docbuild needs quite some RAM (as of May 2018). It sometimes calls
# os.fork() to spawn an external program which then exceeds easily the
# overcommit limit of the system (no RAM is actually used, but this limit is
# very low because there is not even swap on most CI systems.)
ARG MAKEFLAGS_DOCBUILD=$MAKEFLAGS
ENV MAKEFLAGS_DOCBUILD $MAKEFLAGS_DOCBUILD
ARG SAGE_NUM_THREADS_DOCBUILD=$SAGE_NUM_THREADS
ENV SAGE_NUM_THREADS $SAGE_NUM_THREADS_DOCBUILD
RUN make

################################################################################
# Image with a full build of sage, ready to release.                           #
################################################################################
FROM make-all AS make-release
RUN make micro_release

################################################################################
# A releasable (relatively small) image which contains a copy of sage without  #
# temporary build artifacts which is set up to start the command line          #
# interface if no parameters are passed in.                                    #
################################################################################
FROM run-time-dependencies AS sagemath
ARG HOME=/home/sage
ARG SAGE_ROOT=/home/sage/sage
COPY --chown=sage:sage --from=make-release $SAGE_ROOT/ $SAGE_ROOT/
# set path in bashrc to start gap, gp, maxima, ... this does not
# affect Sage (see https://github.com/sagemath/sage/pull/38049)
RUN echo "export PATH=$SAGE_ROOT/local/bin:$PATH" >> ~/.bashrc

COPY ./docker/entrypoint.sh /usr/local/bin/sage-entrypoint
WORKDIR $HOME
ENTRYPOINT ["/usr/local/bin/sage-entrypoint"]
EXPOSE 8888
CMD ["sage"]

################################################################################
# Image with a full build of sage and its documentation but everything         #
# stripped that can be quickly rebuild by make.                                #
################################################################################
FROM make-all AS make-fast-rebuild-clean
# Of course, without site-packages/sage, sage does not start but copying/compiling
# them from build/pkgs/sagelib/src/build is very fast.
RUN make fast-rebuild-clean; \
    rm -rf local/lib/python*/site-packages/sage

################################################################################
# Depending on whether we built from source-clean or not, this image is either #
# identical to make-fast-rebuild-clean or contains a "patch" which can be used #
# to upgrade ARTIFACT_BASE to make-fast-rebuild-clean.                         #
################################################################################
FROM make-fast-rebuild-clean AS sagemath-dev-patch
ARG ARTIFACT_BASE=source-clean
ARG SAGE_ROOT=/home/sage/sage
# Build a patch containing of a tar file which contains all the modified files
# and a list of all modified files (added/updated/removed).
RUN if [ x"$ARTIFACT_BASE" != x"source-clean" ]; then \
        mkdir -p $HOME/patch \
        && find . \( -type f -or -type l \) > $HOME/make-fast-rebuild-clean.manifest \
        && cat $HOME/make-fast-rebuild-clean.manifest $HOME/artifact-base.manifest | sort | uniq -u > $HOME/obsolete \
        && find . \( -type f -or -type l \) -cnewer $HOME/artifact-base.manifest > $HOME/modified \
        && tar -cJf $HOME/patch/modified.tar.xz -T $HOME/modified \
        && cat $HOME/obsolete $HOME/modified | xz > $HOME/patch/modified.xz \
        && rm -rf $SAGE_ROOT \
        && mkdir -p $SAGE_ROOT \
        && mv $HOME/patch $SAGE_ROOT/; \
    fi

################################################################################
# A releasable (relatively small, but still huge) image of this build with all #
# the build artifacts intact so developers can make changes and rebuild        #
# quickly                                                                      #
################################################################################
FROM $ARTIFACT_BASE AS sagemath-dev
ARG SAGE_ROOT=/home/sage/sage
# If docker is backed by aufs, then the following command adds the size of
# ARTIFACT_BASE to the image size. As of mid 2018 this is notably the case with
# the docker instances provided by setup_remote_docker on CircleCI. As a
# result, the sagemath-dev images that are "build-from-latest" are twice AS big
# as the ones that are build on GitLab:
# https://github.com/moby/moby/issues/6119#issuecomment-268870519
COPY --chown=sage:sage --from=sagemath-dev-patch $SAGE_ROOT $SAGE_ROOT
ARG ARTIFACT_BASE=source-clean
# Apply the patch from sagemath-dev-patch if we created one.
RUN if [ x"$ARTIFACT_BASE" != x"source-clean" ]; then \
        echo "Applying `du -hs patch/modified.tar.xz` patch" \
        && xzcat patch/modified.xz | xargs rm -rvf \
        && tar -Jxf patch/modified.tar.xz \
        && rm -rf patch; \
    fi
COPY ./docker/entrypoint-dev.sh /usr/local/bin/sage-entrypoint
ENTRYPOINT ["/usr/local/bin/sage-entrypoint"]
CMD ["bash"]
